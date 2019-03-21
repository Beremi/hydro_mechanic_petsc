import sys, petsc4py
petsc4py.init(sys.argv)

from petsc4py import PETSc

nx, ny, nz = 12, 12, 12       # Dimensions for PETSc DA object

# Set up global space
DMDA3D = PETSc.DMDA().create([nx,ny,nz], dof=3, stencil_width=1)

# Set up a global vector
myVector = DMDA3D.createGlobalVec()
myVector.setName('grid3d')
myVector.set(-1.0)

# Save *.h5 file
ViewHDF5 = PETSc.Viewer()     # Init. Viewer
ViewHDF5.createHDF5('grid.h5', mode=PETSc.Viewer.Mode.WRITE,comm= PETSc.COMM_WORLD)
ViewHDF5.view(obj=myVector)   # Put PETSc object into the viewer
ViewHDF5.destroy()            # Destroy Viewer

# Load *.h5 file
ViewHDF5 = PETSc.Viewer()     # Init. Viewer
ViewHDF5.createHDF5('grid.h5', mode=PETSc.Viewer.Mode.READ,comm= PETSc.COMM_WORLD)
ViewHDF5.view(obj=myVector)   # Put PETSc object into the viewer
ViewHDF5.destroy()            # Destroy Viewer


myVector.destroy()            # Destroy PETSc Vector
DMDA3D.destroy()              # Destroy PETSc DMDA