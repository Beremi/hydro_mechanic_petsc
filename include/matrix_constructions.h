#ifndef MATRIX_CONSTRUCTION_H
#define MATRIX_CONSTRUCTION_H

PetscErrorCode FullMatrix(PetscReal* resMat, const PetscReal h, const PetscReal k);
PetscErrorCode HDivMatrix(PetscReal* resMat, const PetscReal h, const PetscReal k, const PetscReal r);
PetscErrorCode MassMatrix(PetscReal* resMat, const PetscReal h, const PetscReal k);
PetscErrorCode CouplingMatrix(PetscReal* resMat, const PetscReal h);
PetscErrorCode PressureMatrix(PetscReal* resMat, const PetscReal h, const PetscReal pressure);

#endif

