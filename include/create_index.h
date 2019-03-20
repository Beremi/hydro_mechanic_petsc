#ifndef CREATE_INDEX_H
#define CREATE_INDEX_H

PetscErrorCode FullMatrixIndex(const PetscInt i, const PetscInt n, PetscInt* index);
PetscErrorCode MassMatrixIndex(const PetscInt i, const PetscInt n, PetscInt* index);
PetscErrorCode PressureMatrixIndex(const PetscInt i, const PetscInt n, PetscInt* index);
PetscErrorCode FluxBoundaryCondition(PetscInt* p, const PetscInt n, const PetscBool top, const PetscBool bottom, const PetscBool left, const PetscBool right);
PetscErrorCode PressureBoundaryConditionValues(PetscReal* values, const PetscInt n,
		const PetscBool top, const PetscBool bottom, const PetscBool left, const PetscBool right);
#endif

