program mat1
#include <petsc/finclude/petscmat.h>
  use petscmat
  use iso_fortran_env

  implicit none

  PetscErrorCode :: ierr
  PetscInt, parameter :: n = 10
  PetscInt :: i, istart, iend, local_size
  Mat :: A
  PetscScalar, parameter :: one = 1
  Vec :: u, v
  PetscScalar, dimension(:), pointer :: values

  PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER, ierr))

  ! Create identity matrix
  ! - Using MatSetValues:
  !PetscCallA(MatCreate(PETSC_COMM_WORLD, A, ierr))
  !PetscCallA(MatSetType(A, MATMPIAIJ, ierr))
  !PetscCallA(MatSetSizes(A, n, n, PETSC_DETERMINE, PETSC_DETERMINE, ierr))
  !PetscCallA(MatMPIAIJSetPreallocation(A, 1, PETSC_NULL_INTEGER, 0, PETSC_NULL_INTEGER, ierr))
  !
  !PetscCallA(MatGetOwnershipRange(A, istart, iend, ierr))
  !
  !! Using MatDiagonalSet is possible but would require to create a vector full of ones
  !do i = istart,iend-1
  !  PetscCallA(MatSetValues(A, 1, (/i/), 1, (/i/), (/one/), INSERT_VALUES, ierr))
  !enddo
  !
  !PetscCallA(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr))
  !PetscCallA(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr))
  !
  ! - Using MatCreateConstantDiagonal:
  PetscCallA(MatCreateConstantDiagonal(PETSC_COMM_WORLD, n, n, PETSC_DETERMINE, PETSC_DETERMINE, one, A,ierr));

  PetscCallA(MatScale(A, real(0.1, kind=kind(one)), ierr))

  PetscCallA(MatCreateVecs(A, u, v, ierr))

  PetscCallA(VecGetOwnershipRange(u, istart, iend, ierr))
  local_size = iend - istart

  PetscCallA(VecGetArrayF90(u, values, ierr))
  do i = 1,local_size
    values(i) = PetscIntToReal(i + istart) * 10
  enddo
  PetscCallA(VecRestoreArrayF90(u, values, ierr))

  PetscCallA(MatMult(A, u, v, ierr))

  PetscCallA(VecView(v, PETSC_VIEWER_STDOUT_WORLD, ierr))

  PetscCallA(VecDestroy(u, ierr))
  PetscCallA(VecDestroy(v, ierr))
  PetscCallA(MatDestroy(A, ierr))

  PetscCallA(PetscFinalize(ierr))
end program mat1
