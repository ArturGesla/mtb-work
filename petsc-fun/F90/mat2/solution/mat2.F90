program mat2
#include <petsc/finclude/petscmat.h>
  use petscmat
  use iso_fortran_env

  implicit none

  PetscErrorCode :: ierr
  PetscInt, parameter :: NpointsPerDir = 10, N = NpointsPerDir * NpointsPerDir, stencil_size = 5
  PetscInt :: i, istart, iend, colIndex, rowIndex
  PetscScalar, parameter ::  h = 1.0 / (NpointsPerDir + 1)
  Mat :: A
  PetscInt, dimension(stencil_size) :: indexes
  PetscScalar, dimension(stencil_size) :: values = (/ -1.0 / (h * h), -1.0 / (h * h), 4.0 / (h * h), -1.0 / (h * h), -1.0 / (h * h) /)

  PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER, ierr))

  PetscCallA(MatCreate(PETSC_COMM_WORLD, A, ierr))
  PetscCallA(MatSetType(A, MATMPIAIJ, ierr))
  PetscCallA(MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, N, N, ierr))
  ! Comment the next line and uncomment the two lines which follow to check
  ! the effect preallocation has on performance
  PetscCallA(MatMPIAIJSetPreallocation(A, stencil_size, PETSC_NULL_INTEGER, stencil_size - 1, PETSC_NULL_INTEGER, ierr))
  ! PetscCallA(MatMPIAIJSetPreallocation(A, 0, PETSC_NULL_INTEGER, 0 , PETSC_NULL_INTEGER, ierr))
  ! PetscCallA(MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr))

  PetscCallA(MatGetOwnershipRange(A, istart, iend, ierr))

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

  ! Uncomment the next line if you want to view the matrix content
  ! PetscCallA(MatView(A, PETSC_VIEWER_STDOUT_WORLD, ierr))
  ! PetscCallA(MatView(A, PETSC_VIEWER_DRAW_WORLD, ierr)) ! use it with -draw_pause <sec>

  PetscCallA(MatDestroy(A, ierr))

  PetscCallA(PetscFinalize(ierr))
end program mat2
