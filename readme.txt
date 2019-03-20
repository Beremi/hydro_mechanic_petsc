Compilation:
make

Depends on PETSc and correctly set variables for the PETSc makefile. Random permeability generation depends on GSL library.

Typical usage
mpirun -n 4 ./hdiv -mrPrecond 1 -fieldsplit_D_pc_type none -fieldsplit_D_ksp_type gmres -ksp_type gmres -ksp_monitor -ksp_converged_reason -fieldsplit_M_ksp_monitor -fieldsplit_M_ksp_converged_reason -randPerm 1 -k_scale 1e-6 k_sigma 1 -cpp 1e-4 -n 50 -tau 0.1

-mrPrecond Turns on the preconditioning for hdiv block
-asm Sets the number of overlapping layers for Schwarz preconditioner (Note that this means that the length of the layer is asm*1/n)
-randPerm Turns on random generation of permeability
-k_scale Sets the scale of the permeability
-k_sigma Sets sigma parameter for the lognormal distribution used for generation of permeability
-cpp Sets the c_{pp} parameter, storativity
-n Sets the space discretization parameter
-tau Sets the length of the timestep

The mpi parameter n also sets the number of subdomains used for Schwarz method.
