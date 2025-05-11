#include <petscvec.h>
#include <petsctime.h>

int main(int argc, char **argv)
{
  // Use a smaller vector size if the program runs out of memory
  const PetscInt N = 600000000;
  Vec bigVec;
  PetscInt i, istart, iend, local_size;
  PetscInt *indexes = NULL;
  PetscScalar *values = NULL;
  PetscLogDouble t1, t2;

  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));

  PetscCall(VecCreate(PETSC_COMM_WORLD, &bigVec));
  PetscCall(VecSetType(bigVec, VECMPI));
  PetscCall(VecSetSizes(bigVec, PETSC_DECIDE, N));

  PetscCall(VecGetOwnershipRange(bigVec, &istart, &iend));

  PetscCall(PetscTime(&t1));
  local_size = iend - istart;
  indexes = (PetscInt*)malloc(local_size * sizeof(PetscInt));
  values = (PetscScalar*)malloc(local_size * sizeof(PetscScalar));
  for (i = 0; i < local_size; i++) {
    indexes[i] = i + istart;
    values[i] = (PetscScalar)indexes[i];
  }
  PetscCall(VecSetValues(bigVec, local_size, indexes, values, INSERT_VALUES));
  free(indexes);
  free(values);

  PetscCall(VecAssemblyBegin(bigVec));
  PetscCall(VecAssemblyEnd(bigVec));
  PetscCall(PetscTime(&t2));

  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Vector initialized in %lf s\n", t2 - t1));

  PetscCall(VecDestroy(&bigVec));

  PetscCall(PetscFinalize());

  return 0;
}
