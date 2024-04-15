import meshio
import matplotlib.pyplot as plt
from dolfin import * 
#msh = meshio.read("import_stl.msh")
msh=meshio.read("pipe_cad.msh")
meshio.write("pipe_cad.xdmf", meshio.Mesh(points=msh.points, cells={"tetra": msh.cells["tetra"]}))
#boundary face
meshio.write("mf_pipe_cad.xdmf", meshio.Mesh(points=msh.points, cells={"triangle": msh.cells["triangle"]},cell_data={"triangle": {"name_to_read": msh.cell_data["triangle"]["gmsh:physical"]}}))
#boundary cell
meshio.write("cf_pipe_cad.xdmf", meshio.Mesh(points=msh.points, cells={"tetra": msh.cells["tetra"]},cell_data={"tetra": {"name_to_read":msh.cell_data["tetra"]["gmsh:physical"]}}))
mesh = Mesh()
with XDMFFile("pipe_cad.xdmf") as infile:
    infile.read(mesh)
##
mvc = MeshValueCollection("size_t", mesh, 2) 
with XDMFFile("mf_pipe_cad.xdmf") as infile:
    infile.read(mvc, "name_to_read")
mf = cpp.mesh.MeshFunctionSizet(mesh, mvc)
boundary_markers=MeshFunction("size_t",mesh,mvc)
##
mvc = MeshValueCollection("size_t", mesh, 3)
with XDMFFile("cf_pipe_cad.xdmf") as infile:
    infile.read(mvc, "name_to_read")
cf = cpp.mesh.MeshFunctionSizet(mesh, mvc)
T=1
in_marker=1
out_marker=2
wall_marker=3
domain_marker=4
#deltat=T/num_steps
dt=0.01
num_steps=100
mu=0.0035
rho=1060
nu=mu/rho
V=VectorFunctionSpace(mesh,"P",3) #velocity space
Q=FunctionSpace(mesh,"P",1) #pressure space
bcu_in=DirichletBC(V,Constant((0.037386,-0.002018,0.196464)),boundary_markers,in_marker)

bcp_out=DirichletBC(Q,Constant(666.61),boundary_markers,out_marker)
bcu_wall=DirichletBC(V,Constant((0.0,0.0,0.0)),boundary_markers,wall_marker)
bcu=[bcu_in,bcu_wall]
bcp=[bcp_out]
u = TrialFunction(V)
v = TestFunction(V)
p = TrialFunction(Q)
q = TestFunction(Q)

# t=n (u,p) unknown
#u_ = Function(V)
#p_ = Function(Q)

# t=n-1 known
#u_1 = Function(V)
#p_1 = Function(Q)
# Create functions
u0 = Function(V)
u1 = Function(V)
p1 = Function(Q)

# Define coefficients
k = Constant(dt)
f = Constant((0, 0, 0))

# Tentative velocity step
F1 = (1/k)*inner(u - u0, v)*dx +inner(grad(u0)*u0, v)*dx + \
     nu*inner(grad(u), grad(v))*dx - inner(f, v)*dx
a1 = lhs(F1)
L1 = rhs(F1)

# Pressure update
a2 = inner(grad(p), grad(q))*dx
L2 = -(rho/k)*div(u1)*q*dx

# Velocity update
a3 = inner(u, v)*dx
L3 = inner(u1, v)*dx - (k/rho)*inner(grad(p1), v)*dx

# Assemble matrices
A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)

# Use amg preconditioner if available
prec = "amg" if has_krylov_solver_preconditioner("amg") else "default"

# Use nonzero guesses - essential for CG with non-symmetric BC
parameters['krylov_solver']['nonzero_initial_guess'] = True

# Create files for storing solution
ufile = File("results/velocity_cad.pvd")
pfile = File("results/pressure_cad.pvd")

# Time-stepping
t = dt
while t < T + DOLFIN_EPS:
    print(t)
    # Update pressure boundary condition
    #p_in.t = t
    begin("Computing tentative velocity")
    # Compute tentative velocity step
    b1 = assemble(L1)
    [bc.apply(A1, b1) for bc in bcu]
    solve(A1, u1.vector(), b1, "bicgstab", "default")
    end()
    begin("Computing pressure correction")
    # Pressure correction
    b2 = assemble(L2)
    [bc.apply(A2, b2) for bc in bcp]
    [bc.apply(p1.vector()) for bc in bcp]
    solve(A2, p1.vector(), b2, "bicgstab", prec)
    end()
    begin("Computing velocity correction")
    # Velocity correction
    b3 = assemble(L3)
    [bc.apply(A3, b3) for bc in bcu]
    solve(A3, u1.vector(), b3, "bicgstab", "default")
    end()
    # Save to file
    ufile << u1
    pfile << p1

    # Move to next time step
    u0.assign(u1)
    t += dt

# Plot solution
plt.figure()
plot(p1, title="Pressure")

#plt.figure()
#plot(u1, title="Velocity")

plt.show()
