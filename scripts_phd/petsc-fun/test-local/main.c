#include <petsc.h>
int main(int argc, char **argv)
{
printf("hello\n");


PetscMPIInt rank, size;
PetscFunctionBeginUser;
PetscInitialize(&argc,&argv,NULL,NULL);
MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
MPI_Comm_size(PETSC_COMM_WORLD, &size);


//create
PetscInt rankp1 = rank +1;
PetscScalar sizes2 = size*1.0/2.0;
Vec x;
VecCreateMPI(PETSC_COMM_WORLD,rankp1,PETSC_DETERMINE,&x);
VecSet(x,sizes2);
VecView(x,PETSC_VIEWER_STDOUT_WORLD);

/*
//duplicate
Vec y;
(VecDuplicate(x,&y));
(VecCopy(x,y));
*/
/*
PetscScalar *array;
PetscScalar ix[rankp1];
(VecGetArray(x,&array));

for (int i=0; i<rankp1; i++)
{
ix[i]=i;

printf("value is: %d\t%lf\n",ix[i],array[i]);
}
(VecSetValues(y,rankp1,&ix,&array,INSERT_VALUES));

*/
/*
(VecView(y,PETSC_VIEWER_STDOUT_WORLD));

//PetscScalar val;
//(VecDot(x,y,&val));
//printf("dot prod:%lf",val);


PetscScalar product,norm;

(VecDot(x,y,&product));
(VecNorm(y,NORM_2,&norm));
(PetscPrintf(PETSC_COMM_WORLD,"Hello: %lf\t%lf\n",product,norm*norm));

*/
//(PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Rank %d out of %d says hello \n", rank, size ));
PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);

//(VecDestroy())
PetscFinalize();

  return 0;
}
