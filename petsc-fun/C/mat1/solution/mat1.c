#include <petscmat.h>

int main(int argc, char **argv)
{
  const PetscInt n = 10;
  PetscInt i, istart, iend, local_size;
  Mat A;
  const PetscScalar one = 1.0;
  Vec u, v;
  PetscScalar *values = NULL;

  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));

  // Create identity matrix
  /* - Using MatSetValues:
  PetscCall(MatCreate(PETSC_COMM_WORLD, &A));
  PetscCall(MatSetType(A, MATMPIAIJ));
  PetscCall(MatSetSizes(A, n, n, PETSC_DETERMINE, PETSC_DETERMINE));
  PetscCall(MatMPIAIJSetPreallocation(A, 1, NULL, 0, NULL));

  PetscCall(MatGetOwnershipRange(A, &istart, &iend));

  // Using MatDiagonalSet is possible but would require to create a vector full of ones
  for (i = istart; i < iend; i++) {
    PetscCall(MatSetValues(A, 1, &i, 1, &i, &one, INSERT_VALUES));
  }

  PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
  */
  // - Using MatCreateConstantDiagonal:
  PetscCall(MatCreateConstantDiagonal(PETSC_COMM_WORLD, n, n, PETSC_DETERMINE, PETSC_DETERMINE, 1., &A));

  PetscCall(MatScale(A, 0.1));

  PetscCall(MatCreateVecs(A, &u, &v));

  PetscCall(VecGetOwnershipRange(u, &istart, &iend));
  local_size = iend - istart;

  PetscCall(VecGetArray(u, &values));
  for (i = 0; i < local_size; i++) {
    values[i] = (PetscScalar)(i + istart + 1) * 10;
  }
  PetscCall(VecRestoreArray(u, &values));

  PetscCall(MatMult(A, u, v));

  PetscCall(VecView(v, PETSC_VIEWER_STDOUT_WORLD));

  PetscCall(VecDestroy(&u));
  PetscCall(VecDestroy(&v));
  PetscCall(MatDestroy(&A));

  PetscCall(PetscFinalize());

  return 0;
}
