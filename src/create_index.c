#include <petscsys.h>

PetscErrorCode FullMatrixIndex(const PetscInt i, const PetscInt n, PetscInt* index){
  PetscInt start;

  PetscFunctionBeginUser;
  // x
  start = i + (n+1)*(i/n);
  index[0] = start;
  index[1] = start + 2*n + 1;
  // y
  index[2] = start + n;
  index[3] = start + n + 1;
  // p
  index[4] = 2*n*(n+1) + i;
  PetscFunctionReturn(0);
}

PetscErrorCode MassMatrixIndex(const PetscInt i, const PetscInt n, PetscInt* index){
  PetscInt start;

  PetscFunctionBeginUser;
  // x
  start = i + (n+1)*(i/n);
  index[0] = start;
  index[1] = start + 2*n + 1;
  // y
  index[2] = start + n;
  index[3] = start + n + 1;
  PetscFunctionReturn(0);
}

PetscErrorCode PressureMatrixIndex(const PetscInt i, const PetscInt n, PetscInt* index){
  PetscFunctionBeginUser;
  // p
  index[0] = i;
  PetscFunctionReturn(0);
}

PetscErrorCode FluxBoundaryCondition(PetscInt* index, const PetscInt n,
		const PetscBool top, const PetscBool bottom, const PetscBool left, const PetscBool right){
  PetscInt i;
  PetscInt shift = 0;

  PetscFunctionBeginUser;
  if(top){
	  for(i = 0; i < n; i++){
		    index[i+shift] = 2*n*n + n + i;
	  }
	  shift += n;
  }
  if(bottom){
	  for(i = 0; i < n; i++){
		    index[i+shift] = i;
	  }
	  shift += n;
  }
  if(left){
	  for(i = 0; i < n; i++){
		    index[i+shift] = n + i*(2*n + 1);
	  }
	  shift += n;
  }
  if(right){
	  for(i = 0; i < n; i++){
		    index[i+shift] = 2*n + i*(2*n + 1);
	  }
	  shift += n;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode PressureBoundaryConditionValues(PetscReal* values, const PetscInt n,
		const PetscBool top, const PetscBool bottom, const PetscBool left, const PetscBool right){
  PetscInt i;
  PetscInt shift = 0;
  PetscReal h = 1.0/n;

  PetscFunctionBeginUser;
  if(top){
	  for(i = 0; i < n; i++){
		  values[i+shift] = h;
	  }
	  shift += n;
  }
  if(bottom){
	  for(i = 0; i < n; i++){
		  values[i+shift] = 0;
	  }
	  shift += n;
  }
  if(left){
	  for(i = 0; i < n; i++){
		  values[i+shift] = -(i*h*h+h*h/2);
	  }
	  shift += n;
  }
  if(right){
	  for(i = 0; i < n; i++){
		  values[i+shift] = i*h*h+h*h/2;
	  }
	  shift += n;
  }
  PetscFunctionReturn(0);
}

