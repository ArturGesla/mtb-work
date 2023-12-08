program vec2a
#include <petsc/finclude/petscvec.h>
  use petscvec
  use iso_c_binding
  use iso_fortran_env

  implicit none

  ! Use a smaller vector size if the program runs out of memory
  PetscInt, parameter :: N = 600000000
  PetscErrorCode :: ierr
  Vec :: bigVec
  PetscInt :: i, istart, iend
  PetscLogDouble :: t1, t2
  character(len=80) :: str

  PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER, ierr))

  PetscCallA(VecCreate(PETSC_COMM_WORLD, bigVec, ierr))
  PetscCallA(VecSetType(bigVec, VECMPI, ierr))
  PetscCallA(VecSetSizes(bigVec, PETSC_DECIDE, N, ierr))

  PetscCallA(VecGetOwnershipRange(bigVec, istart, iend, ierr))

  PetscCallA(PetscTime(t1, ierr))
  do i = istart,iend-1
    PetscCallA(VecSetValue(bigVec, i, PetscIntToReal(i), INSERT_VALUES, ierr))
  enddo

  PetscCallA(VecAssemblyBegin(bigVec, ierr))
  PetscCallA(VecAssemblyEnd(bigVec, ierr))
  PetscCallA(PetscTime(t2, ierr))

  write(str, "('Vector initialized in ',F0.6,' s',A)") t2 - t1, C_NEW_LINE
  PetscCallA(PetscPrintf(PETSC_COMM_WORLD, trim(str), ierr))

  PetscCallA(VecDestroy(bigVec, ierr))

  PetscCallA(PetscFinalize(ierr))
end program vec2a
