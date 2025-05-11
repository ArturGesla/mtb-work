#include <petscmat.h>

int main(int argc, char **argv)
{
  const PetscInt NpointsPerDir = 10, N = NpointsPerDir * NpointsPerDir, stencil_size = 5;
  PetscInt i, istart, iend, colIndex, rowIndex;
  const PetscScalar h = 1.0 / (NpointsPerDir + 1);
  Mat A;
  PetscInt indexes[stencil_size];
  PetscScalar values[] = { -1.0 / (h * h), -1.0 / (h * h), 4.0 / (h * h), -1.0 / (h * h), -1.0 / (h * h) };

  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));
  PetscCall(MatCreate(PETSC_COMM_WORLD, &A));
  PetscCall(MatSetType(A, MATMPIAIJ));
  PetscCall(MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, N, N));
  // Comment the next line and uncomment the two lines which follow to check
  // the effect preallocation has on performance
  PetscCall(MatMPIAIJSetPreallocation(A, stencil_size, NULL, stencil_size - 1, NULL));
  // PetscCall(MatMPIAIJSetPreallocation(A, 0, NULL, 0 , NULL));
  // PetscCall(MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));

  PetscCall(MatGetOwnershipRange(A, &istart, &iend));

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

  // Uncomment the next line if you want to view the matrix content
  // PetscCall(MatView(A, PETSC_VIEWER_STDOUT_WORLD));
  // PetscCall(MatView(A, PETSC_VIEWER_DRAW_WORLD)); // use it with -draw_pause <sec>

  PetscCall(MatDestroy(&A));

  PetscCall(PetscFinalize());

  return 0;
}
