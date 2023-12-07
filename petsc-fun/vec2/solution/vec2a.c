#include <petscvec.h>
#include <petsctime.h>

int main(int argc, char **argv)
{
  // Use a smaller vector size if the program runs out of memory
  const PetscInt N = 600000000;
  Vec bigVec;
  PetscInt i, istart, iend;
  PetscLogDouble t1, t2;

  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));

  PetscCall(VecCreate(PETSC_COMM_WORLD, &bigVec));
  PetscCall(VecSetType(bigVec, VECMPI));
  PetscCall(VecSetSizes(bigVec, PETSC_DECIDE, N));

  PetscCall(VecGetOwnershipRange(bigVec, &istart, &iend));

  PetscCall(PetscTime(&t1));
  for (i = istart; i < iend; i++) {
    PetscCall(VecSetValue(bigVec, i, (PetscScalar)i, INSERT_VALUES));
  }

  PetscCall(VecAssemblyBegin(bigVec));
  PetscCall(VecAssemblyEnd(bigVec));
  PetscCall(PetscTime(&t2));

  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Vector initialized in %lf s\n", t2 - t1));

  PetscCall(VecDestroy(&bigVec));

  PetscCall(PetscFinalize());

  return 0;
}
