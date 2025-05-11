#include <petscksp.h>
#include <petsctime.h>

int main(int argc, char **argv)
{
  const PetscInt stencil_size = 5;
  PetscInt NpointsPerDir = 4, N;
  PetscInt i, istart, iend, its, colIndex, rowIndex;
  PetscInt indexes[stencil_size];
  Mat A;
  PetscScalar h, values[stencil_size];
  Vec x, xExact, b, r, error;
  KSP solver;
  const PetscReal rtol = 1.e-8;
  KSPConvergedReason reason;
  PetscReal errorNorm, xnorm, rnorm, bnorm;
  PetscLogDouble t1, t2;
  PetscLogStage stage1, stage2;

  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-size", &NpointsPerDir, NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "NpointsPerDir: %d\n", NpointsPerDir));

  N = NpointsPerDir * NpointsPerDir;
  h = 1.0 / (NpointsPerDir + 1);
  values[0] = values[1] = values[3] = values[4] = -1.0 / (h * h);
  values[2] = 4.0 / (h * h);

  PetscCall(PetscLogStageRegister("Fill & assemble", &stage1));
  PetscCall(PetscLogStageRegister("Solve", &stage2));

  PetscCall(MatCreate(PETSC_COMM_WORLD, &A));
  PetscCall(MatSetType(A, MATMPIAIJ));
  PetscCall(MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, N, N));
  PetscCall(MatMPIAIJSetPreallocation(A, stencil_size, NULL, stencil_size - 1, NULL));

  PetscCall(MatGetOwnershipRange(A, &istart, &iend));

  PetscCall(PetscLogStagePush(stage1));
  for (i = istart; i < iend; i++) {
    rowIndex = i / NpointsPerDir; colIndex = i - rowIndex * NpointsPerDir;

    // Negative indexes are ignored by MatSetValues
    indexes[0] = (rowIndex > 0) ? (i - NpointsPerDir) : -1;
    indexes[1] = (colIndex > 0) ? (i - 1) : -1;
    indexes[2] = i;
    indexes[3] = (colIndex < NpointsPerDir - 1) ? (i + 1) : -1;
    indexes[4] = (rowIndex < NpointsPerDir - 1) ? (i + NpointsPerDir) : -1;

    PetscCall(MatSetValues(A, 1, &i, stencil_size, indexes, values, INSERT_VALUES));
  }

  PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
  PetscCall(PetscLogStagePop());

  // Uncomment the next line (or -mat_view at runtime) if you want to view the matrix content
  // PetscCall(MatView(A, PETSC_VIEWER_STDOUT_WORLD));

  PetscCall(MatCreateVecs(A, &x, &b));
  PetscCall(VecDuplicate(x, &xExact));
  PetscCall(VecSetRandom(xExact, NULL));
  PetscCall(MatMult(A, xExact, b));

  PetscCall(KSPCreate(PETSC_COMM_WORLD, &solver));
  PetscCall(KSPSetOperators(solver, A, A));
  PetscCall(KSPSetTolerances(solver, rtol, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT));

  // Use this (or -ksp_norm_type unpreconditioned at runtime) to compare with the residual computed below.
  // PetscCall(KSPSetNormType(solver, KSP_NORM_UNPRECONDITIONED));

  PetscCall(KSPSetFromOptions(solver));
  PetscCall(KSPSetUp(solver));

  // Use this (or -ksp_view at runtime) to view your solver
  PetscCall(KSPView(solver, PETSC_VIEWER_STDOUT_WORLD));

  PetscCall(PetscLogStagePush(stage2));
  PetscCall(PetscTime(&t1));
  PetscCall(KSPSolve(solver, b, x));
  PetscCall(PetscTime(&t2));
  PetscCall(PetscLogStagePop());

  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Elapse time: %lf s\n", t2 - t1));

  PetscCall(KSPGetConvergedReason(solver, &reason));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Converged reason: %s\n", KSPConvergedReasons[reason]));

  PetscCall(KSPGetIterationNumber(solver, &its));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Iterations number: %d\n", its));

  PetscCall(KSPGetResidualNorm(solver, &rnorm));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Residual norm: %g\n", rnorm));

  PetscCall(VecDuplicate(x, &error));
  PetscCall(VecCopy(x, error));
  PetscCall(VecAXPY(error, -1.0, xExact));
  PetscCall(VecNorm(error, NORM_2, &errorNorm));
  PetscCall(VecNorm(x, NORM_2, &xnorm));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Error norm (NORM_2): %g, relative error: %g\n", errorNorm, errorNorm / xnorm));
  PetscCall(VecNorm(error, NORM_INFINITY, &errorNorm));
  PetscCall(VecNorm(x, NORM_INFINITY, &xnorm));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Error norm (NORM_INFINITY): %g, relative error: %g\n", errorNorm, errorNorm / xnorm));

  PetscCall(VecNorm(b, NORM_2, &bnorm));
  PetscCall(VecDuplicate(b, &r));
  PetscCall(MatMult(A, x, r));
  PetscCall(VecAXPY(r, -1.0, b));
  PetscCall(VecNorm(r, NORM_2, &rnorm));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Residual norm (computed): %g, relative error: %g\n", rnorm, rnorm / bnorm));

  PetscCall(KSPDestroy(&solver));
  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&xExact));
  PetscCall(VecDestroy(&b));
  PetscCall(VecDestroy(&r));
  PetscCall(VecDestroy(&error));
  PetscCall(MatDestroy(&A));

  PetscCall(PetscFinalize());

  return 0;
}
