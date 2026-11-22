##############################################
#  Diffusion coupled with linear elasticity  #
#            a mixed formulation             #
##############################################

# IMPORTANT: This code has not been validated
#            and may contain errors 
#
# for validation, please compare to [M.Ahmadi and B.Eidel (2019)] 

# Import Firedrake
from firedrake import *

# Import matplotlib
import matplotlib.pyplot as plt

## Input parameters
# Material properties
E     = Constant(80e9)     # Young's modulus
nu    = Constant(0.22)     # Poisson's ratio
N_A   = Constant(6.022e23) # Avogadro's number
k_B   = Constant(1.38e-23) # Boltzmann's constant
T     = Constant(298.15)   # Absolute temperature
M     = Constant(500)      # Molecular mobility
Omega = Constant(8.5e-6)   # Partial molar volume
c_max = Constant(88.67e3)  # Maximum concentration
c_0   = Constant(1.0e3)    # Initial concentration
# Time parameters
dt      = 0.025         # Time step size
t_end   = 50.0          # Final time
n_steps = int(t_end/dt) # Number of time steps

# Import mesh from Gmsh
mesh = Mesh("mesh1.msh")
# show the mesh
fig, axes = plt.subplots()
triplot(mesh, axes=axes)
plt.show()

# Elasticity parameters
lmbda   = E*nu / ((1 + nu)*(1 - 2*nu))  # First Lame parameter
mu      = E / (2*(1 + nu))              # Second Lame parameter (shear modulus)

## Strain and Stress tensors
I = Identity(3)  # 3D identity tensor for plane strain

def epsilon(u):
    return sym(grad(u))

def epsilon_c(c):
    return (Omega/3.0) * (c - c_0) * I

def epsilon_el(u, c):
    return epsilon(u) - epsilon_c(c)

def sigma(u, c):
    return 2*mu*epsilon_el(u, c) + lmbda*tr(epsilon_el(u, c))*I

def sigma_2D(u, c):
    sig = sigma(u, c)
    return as_tensor(((sig[0, 0], sig[0, 1]),
                      (sig[1, 0], sig[1, 1])))


## FE settings
# Function space (P2 for displacement, P1 for concentration)
V     = VectorFunctionSpace(mesh, "CG", 2)  # Quadratic for displacement
W     = FunctionSpace(mesh, "CG", 1)        # Linear for concentration
Z     = V * W                               # Mixed function space
V_out = VectorFunctionSpace(mesh, "CG", 1)  # Piecewise linear for output

# Trial and Test functions
z    = Function(Z)      # (u, c) pair
u, c = split(z)         # Extract displacement and concentration fields
z_old = Function(Z)     
_, c_old = split(z_old) # Extract old concentration field 
del_u, del_c = TestFunctions(Z)  # Test functions


## Boundary conditions
# Dirichlet boundary conditions for displacement
# Fix left and bottom in x and y directions
ubar  = Constant(0.0) # Displacement on left and bottom
bcs_u = [DirichletBC(Z.sub(0).sub(0), ubar, 2),  # u_x fixed on left
         DirichletBC(Z.sub(0).sub(1), ubar, 3)]  # u_y fixed on bottom
# Dirichlet boundary conditions for concentration
cbar  = c_max # Concentration on arc
bcs_c = [DirichletBC(Z.sub(1), cbar, 1)]
# Combine boundary conditions
bcs = bcs_u + bcs_c

## Initial conditions
z_old.sub(1).interpolate(c_0)
z.sub(1).interpolate(c_0)

## Weak Forms
# 1. Elasticity
R_u = inner(grad(del_u), sigma_2D(u, c)) * dx
# 2. Diffusion
sigma_p = (1.0/3.0) * tr(sigma(u, c)) 
R_c = (c - c_old)/dt * del_c * dx + M * k_B * T * inner(grad(c), grad(del_c)) * dx - (M * c * Omega / N_A) * inner(grad(sigma_p), grad(del_c)) * dx
# Monolithic residual
R = R_u + R_c

# Create the VTK file for output
outfile = VTKFile("solution.pvd")
outfile.write(project(u, V_out, name="u"), project(c, W, name="c"), time=0.0)

# Print problem information
print("======================================")
print("Diffusion coupled with linear elasticity")
print("Number of elements: ", mesh.num_cells())
print(f"Number of DOFs:", Z.dim())
print("======================================")

## Solve the problem
# Newton-Raphson iteration
problem = NonlinearVariationalProblem(R, z, bcs=bcs)
solver  = NonlinearVariationalSolver(problem, solver_parameters={"snes_type": "newtonls"})

for i in range(1, n_steps+1):

    # Solve the system of equations
    solver.solve()
    
    # Update concentration
    z_old.sub(1).assign(z.sub(1))
    
    # Save results
    outfile.write(project(u, V_out, name="u"), project(c, W, name="c"), time=float(dt)*float(i))

    # Print progress
    print(f"time={float(dt)*float(i):.2f}  completed: {100.0*float(i)/float(n_steps):.1f}%")

# Print completion message
print("======================================")
print("Simulation completed successfully \U00002705")
print("======================================")