#ifndef ASM_SHELL_H
#define ASM_SHELL_H

typedef struct {
  PetscInt n,overlap;
  const PetscInt *idx;
  IS isrow;
  VecScatter scatter;
  Vec x,y;
  KSP subksp;
} ASMPC;


PetscErrorCode ASMPCCreate(ASMPC **shell,PetscInt sizeX);
PetscErrorCode ASMPCSetUP(PC pc);
PetscErrorCode ASMPCApply(PC pc, Vec x, Vec y);
PetscErrorCode ASMPCDestroy(PC pc);

#endif

