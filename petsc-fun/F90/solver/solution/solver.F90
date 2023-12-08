program solve
#include <petsc/finclude/petscksp.h>
  use petscksp
  use iso_c_binding
  use iso_fortran_env

  implicit none

  PetscErrorCode :: ierr
  PetscInt, parameter :: stencil_size = 5
  PetscInt :: NpointsPerDir = 4, N
  PetscInt :: i, istart, iend, its, colIndex, rowIndex
  PetscInt, dimension(stencil_size) :: indexes
  Mat :: A
  PetscScalar :: h
  PetscScalar, dimension(stencil_size) :: values
  Vec :: x, xExact, b, r, error
  KSP :: solver
  PetscReal, parameter :: rtol = 1.e-8
  KSPConvergedReason :: reason
  PetscReal :: errorNorm, xnorm, rnorm, bnorm
  PetscLogDouble :: t1, t2
  PetscLogStage :: stage1, stage2
  PetscBool :: flg
  character(len=120) :: str

  PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER, ierr))

  PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-size", NpointsPerDir, flg, ierr))
  write(str, "('NpointsPerDir: ',I4,A)") NpointsPerDir, C_NEW_LINE
  PetscCallA(PetscPrintf(PETSC_COMM_WORLD, trim(str), ierr))

  N = NpointsPerDir * NpointsPerDir
  h = 1.0 / (NpointsPerDir + 1)
  values(:) = -1.0 / (h * h)
  values(3) = 4.0 / (h * h)

  PetscCallA(PetscLogStageRegister("Fill & assemble", stage1, ierr))
  PetscCallA(PetscLogStageRegister("Solve", stage2, ierr))

  PetscCallA(MatCreate(PETSC_COMM_WORLD, A, ierr))
  PetscCallA(MatSetType(A, MATMPIAIJ, ierr))
  PetscCallA(MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, N, N, ierr))
  PetscCallA(MatMPIAIJSetPreallocation(A, stencil_size, PETSC_NULL_INTEGER, stencil_size - 1, PETSC_NULL_INTEGER, ierr))

  PetscCallA(MatGetOwnershipRange(A, istart, iend, ierr))

  PetscCallA(PetscLogStagePush(stage1, ierr))
  do i = istart,iend-1
    rowIndex = i / NpointsPerDir; colIndex = i - rowIndex * NpointsPerDir

    ! Negative indexes are ignored by MatSetValues
    indexes(1) = merge(i - NpointsPerDir, -1, rowIndex > 0)
    indexes(2) = merge(i - 1, -1, colIndex > 0)
    indexes(3) = i
    indexes(4) = merge(i + 1, -1, colIndex < NpointsPerDir - 1)
    indexes(5) = merge(i + NpointsPerDir, -1, rowIndex < NpointsPerDir - 1)

    PetscCallA(MatSetValues(A, 1, i, stencil_size, indexes, values, INSERT_VALUES, ierr))
  enddo

  PetscCallA(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr))
  PetscCallA(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr))
  PetscCallA(PetscLogStagePop(ierr))

  ! Uncomment the next line (or -mat_view at runtime) if you want to view the matrix content
  ! PetscCallA(MatView(A, PETSC_VIEWER_STDOUT_WORLD, ierr))

  PetscCallA(MatCreateVecs(A, x, b, ierr))
  PetscCallA(VecDuplicate(x, xExact, ierr))
  PetscCallA(VecSetRandom(xExact, PETSC_NULL_RANDOM, ierr))
  PetscCallA(MatMult(A, xExact, b, ierr))

  PetscCallA(KSPCreate(PETSC_COMM_WORLD, solver, ierr))
  PetscCallA(KSPSetOperators(solver, A, A, ierr))
  PetscCallA(KSPSetTolerances(solver, rtol, PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER, ierr))

  ! Use this (or -ksp_norm_type unpreconditioned at runtime) to compare with the residual computed below.
  ! PetscCallA(KSPSetNormType(solver, KSP_NORM_UNPRECONDITIONED, ierr))

  PetscCallA(KSPSetFromOptions(solver, ierr))
  PetscCallA(KSPSetUp(solver, ierr))

  ! Use this (or -ksp_view at runtime) to view your solver
  PetscCallA(KSPView(solver, PETSC_VIEWER_STDOUT_WORLD, ierr))

  PetscCallA(PetscLogStagePush(stage2, ierr))
  PetscCallA(PetscTime(t1, ierr))
  PetscCallA(KSPSolve(solver, b, x, ierr))
  PetscCallA(PetscTime(t2, ierr))
  PetscCallA(PetscLogStagePop(ierr))

  write(str, "('Elapse time: ',F0.6,' s',A)") t2 - t1, C_NEW_LINE
  PetscCallA(PetscPrintf(PETSC_COMM_WORLD, trim(str), ierr))

  PetscCallA(KSPGetConvergedReason(solver, reason, ierr))
  write(str, "('Converged reason: ',I4,A)") reason, C_NEW_LINE
  PetscCallA(PetscPrintf(PETSC_COMM_WORLD, trim(str), ierr))

  PetscCallA(KSPGetIterationNumber(solver, its, ierr))
  write(str, "('Iterations number: ',I4,A)") its, C_NEW_LINE
  PetscCallA(PetscPrintf(PETSC_COMM_WORLD, trim(str), ierr))

  PetscCallA(KSPGetResidualNorm(solver, rnorm, ierr))
  write(str, "('Residual norm: ',G12.5,A)") rnorm, C_NEW_LINE
  PetscCallA(PetscPrintf(PETSC_COMM_WORLD, trim(str), ierr))

  PetscCallA(VecDuplicate(x, error, ierr))
  PetscCallA(VecCopy(x, error, ierr))
  PetscCallA(VecAXPY(error, PetscIntToReal(-1), xExact, ierr))
  PetscCallA(VecNorm(error, NORM_2, errorNorm, ierr))
  PetscCallA(VecNorm(x, NORM_2, xnorm, ierr))
  write(str, "('Error norm (NORM_2): ',G12.5,', relative error: ',G12.5,A)") errorNorm, errorNorm / xnorm, C_NEW_LINE
  PetscCallA(PetscPrintf(PETSC_COMM_WORLD, trim(str), ierr))
  PetscCallA(VecNorm(error, NORM_INFINITY, errorNorm, ierr))
  PetscCallA(VecNorm(x, NORM_INFINITY, xnorm, ierr))
  write(str, "('Error norm (NORM_INFINITY): ',G12.5,', relative error: ',G12.5,A)") errorNorm, errorNorm / xnorm, C_NEW_LINE
  PetscCallA(PetscPrintf(PETSC_COMM_WORLD, trim(str), ierr))

  PetscCallA(VecNorm(b, NORM_2, bnorm, ierr))
  PetscCallA(VecDuplicate(b, r, ierr))
  PetscCallA(MatMult(A, x, r, ierr))
  PetscCallA(VecAXPY(r, PetscIntToReal(-1), b, ierr))
  PetscCallA(VecNorm(r, NORM_2, rnorm, ierr))
  write(str, "('Residual norm (computed): ',G12.5,', relative error: ',G12.5,A)") rnorm, rnorm / bnorm, C_NEW_LINE
  PetscCallA(PetscPrintf(PETSC_COMM_WORLD, trim(str), ierr))

  PetscCallA(KSPDestroy(solver, ierr))
  PetscCallA(VecDestroy(x, ierr))
  PetscCallA(VecDestroy(xExact, ierr))
  PetscCallA(VecDestroy(b, ierr))
  PetscCallA(VecDestroy(r, ierr))
  PetscCallA(VecDestroy(error, ierr))
  PetscCallA(MatDestroy(A, ierr))

  PetscCallA(PetscFinalize(ierr))
end program solve
