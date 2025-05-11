#include <petscvec.h>
#include <petsctime.h>

int main(int argc, char **argv)
{
  // Use a smaller vector size if the program runs out of memory
  const PetscInt N = 600000000;
  Vec bigVec;
  PetscInt i, istart, iend, local_size;
  PetscScalar *values = NULL;
  PetscLogDouble t1, t2;

  PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));

  PetscCall(VecCreate(PETSC_COMM_WORLD, &bigVec));
  PetscCall(VecSetType(bigVec, VECMPI));
  PetscCall(VecSetSizes(bigVec, PETSC_DECIDE, N));

  PetscCall(VecGetOwnershipRange(bigVec, &istart, &iend));

  PetscCall(PetscTime(&t1));
  local_size = iend - istart;
  PetscCall(VecGetArray(bigVec, &values));
  for (i = 0; i < local_size; i++) {
    values[i] = (PetscScalar)(i + istart);
  }
  PetscCall(VecRestoreArray(bigVec, &values));
  PetscCall(PetscTime(&t2));

  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Vector initialized in %lf s\n", t2 - t1));

  PetscCall(VecDestroy(&bigVec));

  PetscCall(PetscFinalize());

  return 0;
}
