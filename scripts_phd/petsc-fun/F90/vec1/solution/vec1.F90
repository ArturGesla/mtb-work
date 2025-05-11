program vec1
#include <petsc/finclude/petscvec.h>
  use petscvec
  use iso_c_binding
  use iso_fortran_env

  implicit none

  PetscErrorCode :: ierr
  integer :: size, rank
  Vec :: vec, vec2
  PetscScalar :: product
  PetscReal :: norm
  character(len=120) :: str

  PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER, ierr))

  PetscCallMPIA(MPI_Comm_size(PETSC_COMM_WORLD, size, ierr))
  PetscCallMPIA(MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr))

  PetscCallA(VecCreate(PETSC_COMM_WORLD, vec, ierr))
  PetscCallA(VecSetType(vec, VECMPI, ierr))
  PetscCallA(VecSetSizes(vec, rank + 1, PETSC_DETERMINE, ierr))
  PetscCallA(VecSet(vec, PetscIntToReal(size) / 2, ierr))

  PetscCallA(VecView(vec, PETSC_VIEWER_STDOUT_WORLD, ierr))

  PetscCallA(VecDuplicate(vec, vec2, ierr))
  PetscCallA(VecCopy(vec, vec2, ierr))

  PetscCallA(VecDot(vec, vec2, product, ierr))
  PetscCallA(VecNorm(vec, NORM_2, norm, ierr))

  write(str, "('[',I0,'] dot product = ',F0.6,' | norm = ',F0.6,' | norm² = ',F0.6,A)") rank, product, norm, norm*norm, C_NEW_LINE
  PetscCallA(PetscPrintf(PETSC_COMM_WORLD, trim(str), ierr))

  PetscCallA(VecDestroy(vec, ierr))
  PetscCallA(VecDestroy(vec2, ierr))

  PetscCallA(PetscFinalize(ierr))
end program vec1
