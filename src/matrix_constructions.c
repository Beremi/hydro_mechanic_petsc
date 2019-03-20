#include <petscsys.h>

PetscErrorCode FullMatrix(PetscReal* resMat, const PetscReal h, const PetscReal k){
  PetscFunctionBeginUser;
  resMat[0] = 1.0/3*PetscPowReal(h,2)/k;
  resMat[1] = 1.0/6*PetscPowReal(h,2)/k;
  resMat[2] = 0;
  resMat[3] = 0;
  resMat[4] = h;
  resMat[5] = 1.0/6*PetscPowReal(h,2)/k;
  resMat[6] = 1.0/3*PetscPowReal(h,2)/k;
  resMat[7] = 0;
  resMat[8] = 0;
  resMat[9] = -h;
  resMat[10] = 0;
  resMat[11] = 0;
  resMat[12] = 1.0/3*PetscPowReal(h,2)/k;
  resMat[13] = 1.0/6*PetscPowReal(h,2)/k;
  resMat[14] = h;
  resMat[15] = 0;
  resMat[16] = 0;
  resMat[17] = 1.0/6*PetscPowReal(h,2)/k;
  resMat[18] = 1.0/3*PetscPowReal(h,2)/k;
  resMat[19] = -h;
  resMat[20] = h;
  resMat[21] = -h;
  resMat[22] = h;
  resMat[23] = -h;
  resMat[24] = 0;
  PetscFunctionReturn(0);
}

PetscErrorCode HDivMatrix(PetscReal* resMat, const PetscReal h, const PetscReal k, const PetscReal r){
  PetscFunctionBeginUser;
  resMat[0] = 1.0/3*PetscPowReal(h,2)/k + r;
  resMat[1] = 1.0/6*PetscPowReal(h,2)/k - r;
  resMat[2] = r;
  resMat[3] = -r;
  resMat[4] = 1.0/6*PetscPowReal(h,2)/k - r;
  resMat[5] = 1.0/3*PetscPowReal(h,2)/k + r;
  resMat[6] = -r;
  resMat[7] = r;
  resMat[8] = r;
  resMat[9] = -r;
  resMat[10] = 1.0/3*PetscPowReal(h,2)/k + r;
  resMat[11] = 1.0/6*PetscPowReal(h,2)/k - r;
  resMat[12] = -r;
  resMat[13] = r;
  resMat[14] = 1.0/6*PetscPowReal(h,2)/k - r;
  resMat[15] = 1.0/3*PetscPowReal(h,2)/k + r;
  PetscFunctionReturn(0);
}

/* Individual matrices mass M, coupling B^T, pressure -D; sys = [M B^T; B -D] */
PetscErrorCode MassMatrix(PetscReal* resMat, const PetscReal h, const PetscReal k){
  PetscFunctionBeginUser;
  resMat[0] = 1.0/3*PetscPowReal(h,2)/k;
  resMat[1] = 1.0/6*PetscPowReal(h,2)/k;
  resMat[2] = 0;
  resMat[3] = 0;
  resMat[4] = 1.0/6*PetscPowReal(h,2)/k;
  resMat[5] = 1.0/3*PetscPowReal(h,2)/k;
  resMat[6] = 0;
  resMat[7] = 0;
  resMat[8] = 0;
  resMat[9] = 0;
  resMat[10] = 1.0/3*PetscPowReal(h,2)/k;
  resMat[11] = 1.0/6*PetscPowReal(h,2)/k;
  resMat[12] = 0;
  resMat[13] = 0;
  resMat[14] = 1.0/6*PetscPowReal(h,2)/k;
  resMat[15] = 1.0/3*PetscPowReal(h,2)/k;
  PetscFunctionReturn(0);
}

PetscErrorCode CouplingMatrix(PetscReal* resMat, const PetscReal h){
  PetscFunctionBeginUser;
  resMat[0] = h;
  resMat[1] = -h;
  resMat[2] = h;
  resMat[3] = -h;
  PetscFunctionReturn(0);
}

PetscErrorCode PressureMatrix(PetscReal* resMat, const PetscReal h, const PetscReal pressure){
  PetscFunctionBeginUser;
  resMat[0] = PetscPowReal(h,2)*pressure;
  PetscFunctionReturn(0);
}

