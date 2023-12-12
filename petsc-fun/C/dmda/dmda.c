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
  TBC

  // View the DMDA object
  TBC

  // Create the A matrix from the DMDA object
  TBC

  // Retrieve local information from the DMDA object
  TBC

  hx2 = 1.0 / ((info.mx - 1) * (info.mx - 1));
  hy2 = 1.0 / ((info.my - 1) * (info.my - 1));

  coef = 1.0;
  coef5[0] = 2.0 / hx2 + 2.0 / hy2;
  coef5[1] = -1.0 / hx2; coef5[2] = -1.0 / hx2;
  coef5[3] = -1.0 / hy2; coef5[4] = -1.0 / hy2;

  // Loop on the grid points
  for ( TBC ) {

    for ( TBC ) {

      if ( TBC ) {
        // Set matrix values to enforce boundary conditions (homogeneous Dirichlet conditions)
        TBC
      } else {
        // Set matrix values fo interior points
	TBC
      }
    }
  }
  PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
  // PetscCall(MatView(A, PETSC_VIEWER_STDOUT_WORLD));

  // Create global vectors b and x from the DMDA object
  TBC

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
