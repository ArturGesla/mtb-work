#include <petsc.h>
int main(int argc, char **argv)
{
printf("hello\n");


PetscMPIInt rank, size;
PetscFunctionBeginUser;
PetscCall(PetscInitialize(&argc,&argv,NULL,NULL));
PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));


//create
PetscInt rankp1 = rank +1;
PetscScalar sizes2 = size*1.0/2.0;
Vec x;
PetscCall(VecCreateMPI(PETSC_COMM_WORLD,rankp1,PETSC_DETERMINE,&x));
PetscCall(VecSet(x,sizes2));
PetscCall(VecView(x,PETSC_VIEWER_STDOUT_WORLD));

//duplicate
Vec y;
PetscCall(VecDuplicate(x,&y));
PetscCall(VecCopy(x,y));

/*
PetscScalar *array;
PetscScalar ix[rankp1];
PetscCall(VecGetArray(x,&array));

for (int i=0; i<rankp1; i++)
{
ix[i]=i;

printf("value is: %d\t%lf\n",ix[i],array[i]);
}
PetscCall(VecSetValues(y,rankp1,&ix,&array,INSERT_VALUES));

*/

PetscCall(VecView(y,PETSC_VIEWER_STDOUT_WORLD));

//PetscScalar val;
//PetscCall(VecDot(x,y,&val));
//printf("dot prod:%lf",val);


PetscScalar product,norm;

PetscCall(VecDot(x,y,&product));
PetscCall(VecNorm(y,NORM_2,&norm));
PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Hello: %lf\t%lf\n",product,norm*norm));


//PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Rank %d out of %d says hello \n", rank, size ));
PetscCall(PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT));

//PetscCall(VecDestroy())
PetscCall(PetscFinalize());

  return 0;
}
