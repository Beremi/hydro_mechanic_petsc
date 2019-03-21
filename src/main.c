static char help[] = "Test of h(div) preconditioner in PETSc";
#include <petscksp.h>
#include "create_index.h"
#include "matrix_constructions.h"
#include "asm_shell.h"
#include "generate_parameters.h"
//comment
int main(int argc, char **args){
  PetscInt n = 30;
  PetscInt timesteps = 10;
  PetscBool time_evolution = PETSC_TRUE;
  PetscInt i,overlap;

  PetscInt dofs_fluxes,dofs_pressures,dofs_total;

  PetscInt ilo,ihi,rowl,rowh;
  IS isLoc,isGlob;

  PetscReal h;
  PetscReal tau = 0.1;
  PetscReal r = 1e6;
  PetscReal* permeabilities;
  PetscReal mu = 0;
  PetscReal sigma = 1;
  PetscReal scale = 1;
  PetscReal cpp = 1;
  PetscReal local_matrix_mass[16];
  PetscReal local_matrix_pressure[1];
  PetscReal local_matrix_coupling[4];
  PetscInt local_index_mass[4];
  PetscInt local_index_pressure[1];

  PetscBool top, bottom, left, right;
  PetscInt* bc_flux_index = NULL;
  PetscReal* bc_flux_values = NULL;
  PetscInt* bc_pressure_index = NULL;
  PetscReal* bc_pressure_values = NULL;
  PetscInt dofs_bc_flux = 0;
  PetscInt dofs_bc_pressure = 0;

  Mat mats[4],nest,precond,M,B,Bt,mD,PmD,P0;
  PetscInt n_loc=PETSC_DECIDE,iend;
  Vec x,b,b0,b1,b_neumann,b_dirichlet;
  IS isg[2];
  KSP ksp,*subksps,*subksp;
  PC pc,subpc;
  PetscBool useASM=PETSC_FALSE,randPerm=PETSC_FALSE,mrPrecond=PETSC_FALSE;
  PetscMPIInt size;
  PetscErrorCode ierr;

  ierr = PetscInitialize(&argc, &args, (char*)0, help);if (ierr) return ierr;
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  if (size == 1) SETERRQ(PETSC_COMM_WORLD,1,"Sequential run is not supported.");
  ierr = PetscOptionsGetInt(NULL,NULL,"-asm",&overlap,&useASM);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-mrPrecond",&mrPrecond,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-timesteps",&timesteps,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-r",&r,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-randPerm",&randPerm,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-k_scale",&scale,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-k_mu",&mu,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-k_sigma",&sigma,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-tau",&tau,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-cpp",&cpp,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-time_evolution",&time_evolution,NULL);CHKERRQ(ierr);

  dofs_fluxes = 2*n*(n+1);
  dofs_pressures = n*n;
  dofs_total = dofs_pressures + dofs_fluxes;
  h = 1.0/n;

  ierr = PetscPrintf(PETSC_COMM_WORLD,
                    "INFO:\n r: %f\n RandPerm: %d\n n: %d\n dofs_fluxes: %d\n dofs_pressures: %d\n dofs_total: %d\n",
                    r,randPerm,n,dofs_fluxes,dofs_pressures,dofs_total);CHKERRQ(ierr);
  ierr = PetscMalloc1(n*n, &permeabilities);CHKERRQ(ierr);

  if (randPerm) {
    #ifdef HAVE_GSL
    GeneratePermeability(permeabilities, n, mu, sigma, scale);
    #else
    SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Program has to be compiled with GSL library");
    #endif
  } else {
    for( i = 0; i<n*n; i++ ){ permeabilities[i] = 1.0; }
  }

  rowh = n;
  rowl = PETSC_DECIDE;
  ierr = PetscSplitOwnership(PETSC_COMM_WORLD,&rowl,&rowh);CHKERRQ(ierr);
  rowh = 0;
  ierr = MPI_Scan(&rowl,&rowh,1,MPIU_INT,MPI_SUM,PETSC_COMM_WORLD);CHKERRQ(ierr);
  rowl = rowh-rowl;
  ilo = rowl*(2*n+1);
  ihi = rowh*(2*n+1);
  if (rowh == n) ihi = dofs_fluxes;

  ierr = MatCreate(PETSC_COMM_WORLD, &M);CHKERRQ(ierr); /* mass = fluxes */
  ierr = MatCreate(PETSC_COMM_WORLD, &mD);CHKERRQ(ierr); /* pressure */
  ierr = MatCreate(PETSC_COMM_WORLD, &Bt);CHKERRQ(ierr); /* coupling */
  ierr = MatSetSizes(M,ihi-ilo,ihi-ilo,dofs_fluxes,dofs_fluxes);CHKERRQ(ierr);
  ierr = MatSetSizes(mD, PETSC_DECIDE, PETSC_DECIDE, dofs_pressures, dofs_pressures);CHKERRQ(ierr);
  ierr = MatSetSizes(Bt,ihi-ilo,PETSC_DECIDE,dofs_fluxes,dofs_pressures);CHKERRQ(ierr);
  ierr = MatSetUp(M);CHKERRQ(ierr);
  ierr = MatSetUp(mD);CHKERRQ(ierr);
  ierr = MatSetUp(Bt);CHKERRQ(ierr);
  /* TODO improve prealloc */
  ierr = MatSeqAIJSetPreallocation(M, 7, NULL);CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(M, 7, NULL,7,NULL);CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(mD, 1, NULL);CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(mD, 1, NULL,0,NULL);CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(Bt, 2, NULL);CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(Bt, 2, NULL,2,NULL);CHKERRQ(ierr);
  if(time_evolution==PETSC_TRUE){
      ierr = MatCreate(PETSC_COMM_WORLD, &PmD);CHKERRQ(ierr);
      ierr = MatSetSizes(PmD, PETSC_DECIDE, PETSC_DECIDE, dofs_pressures, dofs_pressures);CHKERRQ(ierr);
      ierr = MatSetUp(PmD);CHKERRQ(ierr);
      ierr = MatSeqAIJSetPreallocation(PmD, 1, NULL);CHKERRQ(ierr);
      ierr = MatMPIAIJSetPreallocation(PmD, 1, NULL,0,NULL);CHKERRQ(ierr);
  }
  /* TODO improve parallel assebmly (local dofs only) */
  ierr = PetscSplitOwnership(PETSC_COMM_WORLD,&n_loc,&dofs_pressures);CHKERRQ(ierr);
  ierr = MPI_Scan(&n_loc,&iend,1,MPIU_INT,MPI_SUM,PETSC_COMM_WORLD);CHKERRQ(ierr);
  for( i = iend-n_loc; i<iend; i++ ){
    ierr = MassMatrixIndex(i, n, local_index_mass);CHKERRQ(ierr);
    ierr = MassMatrix(local_matrix_mass, h, permeabilities[i]);CHKERRQ(ierr);
    ierr = MatSetValues(M, 4, local_index_mass, 4, local_index_mass, local_matrix_mass, ADD_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  for( i = iend-n_loc; i<iend; i++ ){
    ierr = MassMatrixIndex(i, n, local_index_mass);CHKERRQ(ierr);
    ierr = PressureMatrixIndex(i, n, local_index_pressure);CHKERRQ(ierr);
    ierr = CouplingMatrix(local_matrix_coupling, h);CHKERRQ(ierr);
    ierr = MatSetValues(Bt, 4, local_index_mass, 1, local_index_pressure, local_matrix_coupling, ADD_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(Bt, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  for( i = iend-n_loc; i<iend; i++ ){
    ierr = PressureMatrixIndex(i, n, local_index_pressure);CHKERRQ(ierr);
    ierr = PressureMatrix(local_matrix_pressure,1.0,1.0);CHKERRQ(ierr);
    ierr = MatSetValues(mD, 1, local_index_pressure, 1, local_index_pressure, local_matrix_pressure, ADD_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(mD, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Bt, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(mD, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatCreateTranspose(Bt,&B);CHKERRQ(ierr);
  ierr = MatSetUp(B);CHKERRQ(ierr);
  if(time_evolution==PETSC_TRUE){
      ierr = MatScale(mD, -h*h*cpp/tau);CHKERRQ(ierr);
      ierr = MatCopy(mD, PmD, DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
      ierr = MatScale(PmD, -1.0);CHKERRQ(ierr);
  }else{
      ierr = MatScale(mD, 1.0/r);CHKERRQ(ierr);
  }

  /* Dirichlet boundary condition setup */
  top = PETSC_FALSE;
  bottom = PETSC_FALSE;
  left = PETSC_TRUE;
  right = PETSC_TRUE;

  if(top){dofs_bc_flux += n;}
  if(bottom){dofs_bc_flux += n;}
  if(left){dofs_bc_flux += n;}
  if(right){dofs_bc_flux += n;}
  ierr = PetscMalloc1(dofs_bc_flux, &bc_flux_values);
  ierr = PetscMalloc1(dofs_bc_flux, &bc_flux_index);
  ierr = FluxBoundaryCondition(bc_flux_index, n, top, bottom, left, right);
  for(i = 0; i < dofs_bc_flux; i++){bc_flux_values[i] = 0;}

  if(time_evolution==PETSC_TRUE){
  mats[0]=M;mats[1]=Bt;
  mats[2]=B;mats[3]=mD;
  }else{
  mats[0]=M;mats[1]=Bt;
  mats[2]=B;mats[3]=NULL;
  }
  ierr = MatCreateNest(PETSC_COMM_WORLD,2,NULL,2,NULL,mats,&nest);CHKERRQ(ierr);

  if (mrPrecond) {
    ierr = MatCreate(PETSC_COMM_WORLD, &P0);CHKERRQ(ierr);
    ierr = MatSetSizes(P0, ihi-ilo, ihi-ilo, dofs_fluxes,dofs_fluxes);CHKERRQ(ierr);
    ierr = MatSetUp(P0);CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(P0, 7, NULL);CHKERRQ(ierr);
    ierr = MatMPIAIJSetPreallocation(P0, 7, NULL,7,NULL);CHKERRQ(ierr);
    for( i = iend-n_loc; i<iend; i++ ){
      ierr = MassMatrixIndex(i, n, local_index_mass);CHKERRQ(ierr);
      if(time_evolution==PETSC_TRUE){
          ierr = HDivMatrix(local_matrix_mass,h,permeabilities[i],tau/cpp);CHKERRQ(ierr);
      }else{
          ierr = HDivMatrix(local_matrix_mass,h,permeabilities[i],r*h*h);CHKERRQ(ierr);
      }
      ierr = MatSetValues(P0,4,local_index_mass,4,local_index_mass,local_matrix_mass, ADD_VALUES);CHKERRQ(ierr);
    }
    ierr = MatAssemblyBegin(P0,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(P0,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    ierr = MatZeroRowsColumns(P0, dofs_bc_flux, bc_flux_index, 1.0, NULL, NULL);

    /* TODO reuse nest IS, drop off-diag blocks */
    if(time_evolution==PETSC_TRUE){
    //horni trojuhelnikove predpodmineni, kladny diagonalni blok
    mats[0]=P0;mats[1]=Bt;
    mats[2]=NULL;mats[3]=mD;
    }else{
    mats[0]=P0;mats[1]=NULL;
    mats[2]=NULL;mats[3]=mD;
    }
    ierr = MatCreateNest(PETSC_COMM_WORLD,2,NULL,2,NULL,mats,&precond);CHKERRQ(ierr);
  } else {
    precond = nest;
  }

  //ierr = VecCreate(PETSC_COMM_WORLD, &x);CHKERRQ(ierr);
  //ierr = VecSetSizes(x, PETSC_DECIDE, dofs_total);CHKERRQ(ierr);
  //ierr = VecSetFromOptions(x);CHKERRQ(ierr);
  //ierr = VecDuplicate(x, &b);CHKERRQ(ierr);
  ierr = MatNestSetVecType(nest,VECNEST);CHKERRQ(ierr); /* NOTE: there is a bug in PCFIELDSPLIT with seq VECNEST */
  ierr = MatCreateVecs(nest,&x,&b);CHKERRQ(ierr);

  ierr = VecSet(b, 0.0);CHKERRQ(ierr); /* rhs */
  ierr = VecSet(x, 2.0);CHKERRQ(ierr); /* serves as an initial guess with -ksp_initial_guess_nonzero */

  /* Neumann boundary condition */
  /* TODO: Nice implementation for boundary condition given by functions */
  /* TODO: Separate variables for indication of the boundary conditions */
  top = PETSC_TRUE;
  bottom = PETSC_TRUE;
  left = PETSC_FALSE;
  right = PETSC_FALSE;

  if(top){dofs_bc_pressure += n;}
  if(bottom){dofs_bc_pressure += n;}
  if(left){dofs_bc_pressure += n;}
  if(right){dofs_bc_pressure += n;}
  ierr = PetscMalloc1(dofs_bc_pressure, &bc_pressure_values);CHKERRQ(ierr);
  ierr = PetscMalloc1(dofs_bc_pressure, &bc_pressure_index);CHKERRQ(ierr);
  ierr = FluxBoundaryCondition(bc_pressure_index, n, top, bottom, left, right);CHKERRQ(ierr);
  //ierr = PressureBoundaryConditionValues(bc_pressure_values, n, top, bottom, left, right);CHKERRQ(ierr);
  for( i = 0; i < n; i++){ bc_pressure_values[i] = h;}//top: p = 1
  for( i = n; i < 2*n; i++){ bc_pressure_values[i] = 0;} //bottom: p = 0

  ierr = VecNestGetSubVec(b,0,&b0);CHKERRQ(ierr);
  ierr = VecNestGetSubVec(b,1,&b1);CHKERRQ(ierr);
  ierr = VecDuplicate(b0,&b_neumann);CHKERRQ(ierr);
  ierr = VecDuplicate(b0,&b_dirichlet);CHKERRQ(ierr);
  /* Dirichlet boundary condition */
  ierr = VecSetValues(b_dirichlet, dofs_bc_flux, bc_flux_index, bc_flux_values, INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(b_dirichlet);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(b_dirichlet);CHKERRQ(ierr);
  ierr = VecScale(b_dirichlet, -1.0);CHKERRQ(ierr);
  ierr = MatMultAdd(M, b_dirichlet, b0, b0);CHKERRQ(ierr);
  ierr = MatMultAdd(B, b_dirichlet, b1, b1);CHKERRQ(ierr);
  ierr = VecDestroy(&b_dirichlet);CHKERRQ(ierr);
  /* Neumann boundary condition */
  ierr = VecSetValues(b_neumann, dofs_bc_pressure, bc_pressure_index, bc_pressure_values, INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(b_neumann);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(b_neumann);CHKERRQ(ierr);
  ierr = VecAXPY(b0, 1.0, b_neumann);CHKERRQ(ierr);
  ierr = VecSetValues(b0, dofs_bc_flux, bc_flux_index, bc_flux_values, INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
  ierr = VecDestroy(&b_neumann);CHKERRQ(ierr);
  ierr = MatZeroRowsColumns(M, dofs_bc_flux, bc_flux_index, 1.0, NULL, NULL);CHKERRQ(ierr);
  ierr = MatZeroRows(Bt, dofs_bc_flux, bc_flux_index, 0.0, NULL, NULL);CHKERRQ(ierr);

  ierr = PetscFree(bc_pressure_values);CHKERRQ(ierr);
  ierr = PetscFree(bc_pressure_index);CHKERRQ(ierr);
  ierr = PetscFree(bc_flux_index);CHKERRQ(ierr);
  ierr = PetscFree(bc_flux_values);CHKERRQ(ierr);

  ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);CHKERRQ(ierr);
  ierr = KSPSetType(ksp, KSPMINRES);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp, nest, precond);CHKERRQ(ierr);

  /* construct PCFieldSplit -- options -fieldsplit_[M|D]_*/
  ierr = MatNestGetISs(nest,isg,NULL);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCFIELDSPLIT);CHKERRQ(ierr);
  if(time_evolution==PETSC_TRUE){
      ierr = PCFieldSplitSetType(pc,PC_COMPOSITE_MULTIPLICATIVE);CHKERRQ(ierr);
  }else{
      ierr = PCFieldSplitSetType(pc,PC_COMPOSITE_ADDITIVE);CHKERRQ(ierr);
  }
  ierr = PCFieldSplitSetIS(pc,"M",isg[0]);CHKERRQ(ierr);
  ierr = PCFieldSplitSetIS(pc,"D",isg[1]);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  ierr = KSPSetUp(ksp);CHKERRQ(ierr);

  if (useASM){
    ierr = PCFieldSplitGetSubKSP(pc,NULL,&subksps);CHKERRQ(ierr);
    ierr = KSPGetPC(subksps[0],&subpc);CHKERRQ(ierr);
    ierr = PCSetType(subpc,PCASM);CHKERRQ(ierr);
    ierr = ISCreateStride(PETSC_COMM_WORLD,ihi-ilo,ilo,1,&isLoc);CHKERRQ(ierr);
    ilo -= overlap*(2*n+1);
    if (ilo < 0) ilo = 0;
    ihi += n + overlap*(2*n+1);
    if (ihi>dofs_fluxes) ihi = dofs_fluxes;
    ierr = ISCreateStride(PETSC_COMM_WORLD,ihi-ilo,ilo,1,&isGlob);CHKERRQ(ierr);
    ierr = PCASMSetLocalSubdomains(subpc,1,&isGlob,NULL);CHKERRQ(ierr);
    ierr = PCASMSetType(subpc,PC_ASM_BASIC);CHKERRQ(ierr);
    ierr = KSPSetType(subksps[0],KSPCG);CHKERRQ(ierr);
    ierr = KSPSetFromOptions(subksps[0]);CHKERRQ(ierr);
    ierr = KSPSetUp(subksps[0]);CHKERRQ(ierr);
    ierr = PCASMGetSubKSP(subpc,NULL,NULL,&subksp);CHKERRQ(ierr);
    ierr = KSPGetPC(*subksp,&subpc);CHKERRQ(ierr);
    ierr = PCSetType(subpc,PCLU);CHKERRQ(ierr);
    ierr = PCSetFromOptions(subpc);CHKERRQ(ierr);
    ierr = KSPSetType(*subksp,KSPPREONLY);CHKERRQ(ierr);
  }

  if(time_evolution==PETSC_FALSE){
      ierr = KSPSolve(ksp, b, x); CHKERRQ(ierr);CHKERRQ(ierr);
  }else{
      Vec x1, bi;
      PetscInt count;
      ierr = VecDuplicate(b, &bi);CHKERRQ(ierr);
      for( i = 0; i < timesteps; i++){
          ierr = PetscPrintf(PETSC_COMM_WORLD,"Timestep %i\n",i);CHKERRQ(ierr);
          ierr = VecCopy(b, bi);CHKERRQ(ierr);
          ierr = VecNestGetSubVec(bi,1,&b1);CHKERRQ(ierr);
          ierr = VecNestGetSubVec(x,1,&x1);CHKERRQ(ierr);
          ierr = MatMultAdd(mD, x1, b1, b1);CHKERRQ(ierr);
          ierr = KSPSolve(ksp, bi, x);CHKERRQ(ierr);
          ierr = KSPGetIterationNumber(ksp, &count);CHKERRQ(ierr);
      }
  }

  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = MatDestroy(&M);CHKERRQ(ierr);
  ierr = MatDestroy(&Bt);CHKERRQ(ierr);
  ierr = MatDestroy(&B);CHKERRQ(ierr);
  ierr = MatDestroy(&mD);CHKERRQ(ierr);
  if (time_evolution==PETSC_TRUE){
      ierr = MatDestroy(&PmD);CHKERRQ(ierr);
  }
  ierr = MatDestroy(&nest);CHKERRQ(ierr);
  if (mrPrecond) {
    ierr = MatDestroy(&P0);CHKERRQ(ierr);
    ierr = MatDestroy(&precond);CHKERRQ(ierr);
  }
  ierr = PetscFree(permeabilities);CHKERRQ(ierr);

  ierr = PetscFinalize();
  return ierr;
}

