program solve
#include <petsc/finclude/petscksp.h>
#include <petsc/finclude/petscdmda.h>
  use petscdmda
  use petscksp
  use iso_c_binding
  use iso_fortran_env

  implicit none

  PetscErrorCode :: ierr
  PetscInt, parameter :: nx = 10, ny=10, stencil_size = 5
  PetscInt :: i, j, nb_values, its
  Mat :: A
  Vec :: b, x, error
  KSP :: solver
  PetscReal, parameter :: rtol = 1.e-8
  KSPConvergedReason :: reason
  PetscReal :: errorNorm, rnorm
  PetscLogDouble :: t1, t2
  DM :: dm
  PetscInt, parameter :: stencilWidth = 1
  PetscInt :: mx, my, xs, ys, xm, ym
  MatStencil, dimension(4) :: row
  MatStencil, dimension(4,stencil_size) :: col5
  PetscScalar, dimension(1) :: coef
  PetscScalar, dimension (stencil_size) :: coef5
  PetscScalar :: hx2, hy2
  character(len=120) :: str
  PetscInt, parameter :: ione = 1
  PetscScalar,dimension(:,:), pointer :: bgrid

  PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER, ierr))

  ! Create the 2-D DMDA object
  PetscCallA(DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_STAR, &
                    nx, ny, PETSC_DECIDE, PETSC_DECIDE, 1, stencilWidth, &
                    PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, dm, ierr))

  PetscCallA(DMSetFromOptions(dm, ierr))
  PetscCallA(DMSetUp(dm, ierr))

  ! View the DMDA object
  PetscCallA(DMView(dm, PETSC_VIEWER_STDOUT_WORLD, ierr))

  ! Create the A matrix from the DMDA object
  PetscCallA(DMCreateMatrix(dm, A, ierr))

  ! Retrieve information (mx,my,xs,ys,xm,ym) from the DMDA object
  PetscCallA(DMDAGetInfo(dm, PETSC_NULL_INTEGER, mx, my, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
                   PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
                   PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr))
  PetscCallA(DMDAGetCorners(dm, xs, ys, PETSC_NULL_INTEGER, xm, ym, PETSC_NULL_INTEGER, ierr))

  hx2 = 1.0 / ((mx - 1) * (mx - 1))
  hy2 = 1.0 / ((my - 1) * (my - 1))

  coef = 1.0
  coef5(1) = 2.0 / hx2 + 2.0 / hy2
  coef5(2) = -1.0 / hx2; coef5(3) = -1.0 / hx2
  coef5(4) = -1.0 / hy2; coef5(5) = -1.0 / hy2

  ! Loop on the grid points
  do j = ys, ys + ym - 1
    do i = xs, xs + xm -1
      row(MatStencil_i) = i; row(MatStencil_j) = j; row(MatStencil_c) = 0

      if (i == 0 .or. i == mx - 1 .or. j == 0 .or. j == my - 1) then
        ! Set matrix values to enforce boundary conditions (homogeneous Dirichlet conditions)
        PetscCallA(MatSetValuesStencil(A, ione, row, ione, row, coef, INSERT_VALUES, ierr))
      else
        ! Set matrix values fo interior points
        col5(MatStencil_i,1) = i; col5(MatStencil_j,1) = j; col5(MatStencil_c,1) = 0
        col5(MatStencil_i,2) = i - 1; col5(MatStencil_j,2) = j; col5(MatStencil_c,2) = 0
        col5(MatStencil_i,3) = i + 1; col5(MatStencil_j,3) = j; col5(MatStencil_c,3) = 0
        col5(MatStencil_i,4) = i; col5(MatStencil_j,4) = j - 1; col5(MatStencil_c,4) = 0
        col5(MatStencil_i,5) = i; col5(MatStencil_j,5) = j + 1; col5(MatStencil_c,5) = 0
        PetscCallA(MatSetValuesStencil(A, ione, row, stencil_size, col5, coef5, INSERT_VALUES, ierr))
      endif
    enddo
  enddo

  PetscCallA(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr))
  PetscCallA(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr))
  !The way to init to 1 all the boundary
  !PetscCallA(MatView(A, PETSC_VIEWER_STDOUT_WORLD, ierr))

  ! Create global vectors b and x from the DMDA object
  PetscCallA(DMCreateGlobalVector(dm, b, ierr))
  PetscCallA(DMCreateGlobalVector(dm, x, ierr))

  PetscCallA(VecSet(b, PetscIntToReal(0), ierr))
  !PetscCallA(DMDAVecGetArrayF90(dm,b,bgrid,ierr))
  !if (ys == 0) then
  !  do i=xs, xs+xm-1
  !    bgrid(i,0) = 1.0
  !  end do
  !end if
  !if (ys+ym == my) then
  !  do i=xs, xs+xm-1
  !    bgrid(i,my-1) = 1.0
  !  end do
  !end if
  !if (xs == 0) then
  !  do j=ys, ys+ym-1
  !    bgrid(0,j) = 1.0
  !  end do
  !end if
  !if (xs+xm == mx) then
  !  do j=ys, ys+ym-1
  !    bgrid(mx-1,j) = 1.0
  !  end do
  !end if
  !PetscCallA(DMDAVecRestoreArrayF90(dm,b,bgrid,ierr))
  !PetscCallA(VecView(b,PETSC_VIEWER_STDOUT_WORLD,ierr))
  PetscCallA(VecSetRandom(x, PETSC_NULL_RANDOM, ierr))

  PetscCallA(KSPCreate(PETSC_COMM_WORLD, solver, ierr))
  PetscCallA(KSPSetOperators(solver, A, A, ierr))
  PetscCallA(KSPSetTolerances(solver, rtol, PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER, ierr))
  PetscCallA(KSPSetInitialGuessNonzero(solver, PETSC_TRUE, ierr))
  PetscCallA(KSPSetFromOptions(solver, ierr))
  PetscCallA(KSPSetUp(solver, ierr))

  PetscCallA(KSPView(solver, PETSC_VIEWER_STDOUT_WORLD, ierr))

  PetscCallA(PetscTime(t1, ierr))
  PetscCallA(KSPSolve(solver, b, x, ierr))
  PetscCallA(PetscTime(t2, ierr))

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

  PetscCallA(VecNorm(x, NORM_2, errorNorm, ierr))
  write(str, "('Error norm (NORM_2): ',G12.5,A)") errorNorm, C_NEW_LINE
  PetscCallA(PetscPrintf(PETSC_COMM_WORLD, trim(str), ierr))
  PetscCallA(VecNorm(x, NORM_INFINITY, errorNorm, ierr))
  write(str, "('Error norm (NORM_INFINITY): ',G12.5,A)") errorNorm, C_NEW_LINE
  PetscCallA(PetscPrintf(PETSC_COMM_WORLD, trim(str), ierr))

  PetscCallA(KSPDestroy(solver, ierr))
  PetscCallA(MatDestroy(A, ierr))
  PetscCallA(VecDestroy(b, ierr))
  PetscCallA(VecDestroy(x, ierr))
  PetscCallA(DMDestroy(dm, ierr))

  PetscCallA(PetscFinalize(ierr))
end program solve
