#include <petscdmda.h>
#include <petscksp.h>
#include <petsctime.h>

int main(int argc, char **argv)
{
  const PetscInt nx = 10, ny = 10, stencil_size = 5;
  PetscInt i, j, its;
  Mat A;
  Vec x, b;
  KSP solver;
  const PetscReal rtol = 1.e-8;
  KSPConvergedReason reason;
  PetscReal errorNorm, rnorm;
  PetscLogDouble t1, t2;
  DM dm;
  DMDALocalInfo info;
  const PetscInt stencilWidth = 1;
  MatStencil row, col5[stencil_size];
  PetscScalar hx2, hy2, coef, coef5[stencil_size];
  PetscScalar **bgrid;

  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));
  // Create the 2-D DMDA object
  PetscCall(DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_STAR,
                      nx, ny, PETSC_DECIDE, PETSC_DECIDE, 1, stencilWidth,
                      NULL, NULL, &dm));

  PetscCall(DMSetFromOptions(dm));
  PetscCall(DMSetUp(dm));

  // View the DMDA object
  PetscCall(DMView(dm, PETSC_VIEWER_STDOUT_WORLD));

  // Create the A matrix from the DMDA object
  PetscCall(DMCreateMatrix(dm, &A));

  // Retrieve local information from the DMDA object
  PetscCall(DMDAGetLocalInfo(dm, &info));

  hx2 = 1.0 / ((info.mx - 1) * (info.mx - 1));
  hy2 = 1.0 / ((info.my - 1) * (info.my - 1));

  coef = 1.0;
  coef5[0] = 2.0 / hx2 + 2.0 / hy2;
  coef5[1] = -1.0 / hx2; coef5[2] = -1.0 / hx2;
  coef5[3] = -1.0 / hy2; coef5[4] = -1.0 / hy2;

  // Loop on the grid points
  for (j = info.ys; j < info.ys + info.ym; j++) {

    for (i = info.xs; i < info.xs + info.xm; i++) {
      row.i = i; row.j = j; row.c = 0;

      if (i == 0 || i == (info.mx - 1) || j == 0 || j == (info.my - 1)) {
        // Set matrix values to enforce boundary conditions (homogeneous Dirichlet conditions)
        PetscCall(MatSetValuesStencil(A, 1, &row, 1, &row, &coef, INSERT_VALUES));
      } else {
        // Set matrix values fo interior points
        col5[0].i = i; col5[0].j = j; col5[0].c = 0;
        col5[1].i = i - 1; col5[1].j = j; col5[1].c = 0;
        col5[2].i = i + 1; col5[2].j = j; col5[2].c = 0;
        col5[3].i = i; col5[3].j = j - 1; col5[3].c = 0;
        col5[4].i = i; col5[4].j = j + 1; col5[4].c = 0;
        PetscCall(MatSetValuesStencil(A, 1, &row, stencil_size, col5, coef5, INSERT_VALUES));
      }
    }
  }
  PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
  // PetscCall(MatView(A, PETSC_VIEWER_STDOUT_WORLD));

  // Create global vectors b and x from the DMDA object
  PetscCall(DMCreateGlobalVector(dm, &b));
  PetscCall(DMCreateGlobalVector(dm, &x));

  PetscCall(VecSet(b, 0.0));
  // The way to init to 1 all the boundary
  //PetscCall(DMDAVecGetArray(dm,b,&bgrid));
  //if (info.ys == 0) {
  //  for (i=info.xs; i<info.xs+info.xm;i++) {
  //    bgrid[0][i] = 1.0; } }
  //if (info.ys+info.ym == info.my) {
  //  for (i=info.xs; i<info.xs+info.xm;i++) {
  //    bgrid[info.my-1][i] = 1.0; } }
  //if (info.xs == 0) {
  //  for (j=info.ys; j<info.ys+info.ym; j++) {
  //    bgrid[j][0] = 1.0; } }
  //if (info.xs+info.xm == info.mx) {
  //  for (j=info.ys; j<info.ys+info.ym; j++) {
  //    bgrid[j][info.mx-1] = 1.0; }}
  //PetscCall(DMDAVecRestoreArray(dm,b,&bgrid));
  //PetscCall(VecView(b,PETSC_VIEWER_STDOUT_WORLD));
  PetscCall(VecSetRandom(x, NULL));

  PetscCall(KSPCreate(PETSC_COMM_WORLD, &solver));
  PetscCall(KSPSetOperators(solver, A, A));
  PetscCall(KSPSetTolerances(solver, rtol, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT));
  PetscCall(KSPSetInitialGuessNonzero(solver, PETSC_TRUE));
  PetscCall(KSPSetFromOptions(solver));
  PetscCall(KSPSetUp(solver));

  PetscCall(KSPView(solver, PETSC_VIEWER_STDOUT_WORLD));

  PetscCall(PetscTime(&t1));
  PetscCall(KSPSolve(solver, b, x));
  PetscCall(PetscTime(&t2));

  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Elapse time: %lf s\n", t2 - t1));

  PetscCall(KSPGetConvergedReason(solver, &reason));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Converged reason: %s\n", KSPConvergedReasons[reason]));

  PetscCall(KSPGetIterationNumber(solver, &its));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Iterations number: %d\n", its));

  PetscCall(KSPGetResidualNorm(solver, &rnorm));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Residual norm: %g\n", rnorm));

  PetscCall(VecNorm(x, NORM_2, &errorNorm));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Error norm (NORM_2): %g\n", errorNorm));
  PetscCall(VecNorm(x, NORM_INFINITY, &errorNorm));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Error norm (NORM_INFINITY): %g\n", errorNorm));

  PetscCall(KSPDestroy(&solver));
  PetscCall(MatDestroy(&A));
  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&b));
  PetscCall(DMDestroy(&dm));

  PetscCall(PetscFinalize());

  return 0;
}
