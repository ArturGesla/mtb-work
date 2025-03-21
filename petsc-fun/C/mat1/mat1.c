#include <petsc.h>
int main(int argc, char **argv)
{
printf("hello\n");



PetscMPIInt rank, size;
PetscFunctionBeginUser;
PetscCall(PetscInitialize(&argc,&argv,NULL,NULL));
PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));

Mat A;
PetscCall(MatCreate(PETSC_COMM_WORLD,&A));

PetscInt n=15;
PetscCall(MatSetType(A,MATMPIAIJ));
PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n));

PetscInt istart,iend;
PetscCall(MatGetOwnershipRange(A,&istart,&iend));
PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD,"s %d e %d  \n", istart, iend ));

PetscInt nloc=iend-istart;
//PetscInt xi[nloc];
//PetscScalar vals[nloc];
PetscScalar one=1;
for (int i=istart; i<iend; i++)
{
//xi[i]=istart+i;
//vals[i]=1;
//printf("rank:%d\t el: %i\n",rank,i);
PetscCall(MatSetValues(A,1,&i,1,&i,&one,INSERT_VALUES));
}

//PetscCall(MatSetValues(A,nloc,&xi,nloc,&xi,&vals,INSERT_VALUES));


PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));

PetscCall(MatSetUp(A));

PetscCall(MatScale(A,20.0));

Vec right,left;
PetscCall(MatCreateVecs(A,&right,&left));

for(int i=istart; i<iend; i++)
PetscCall(VecSetValue(right,i,(i+1)*10.0,INSERT_VALUES));

PetscCall(VecAssemblyBegin(right));
PetscCall(VecAssemblyEnd(right));


PetscCall(MatMult(A,right,left));


PetscCall(MatView(A,PETSC_VIEWER_STDOUT_WORLD));
PetscCall(VecView(right,PETSC_VIEWER_STDOUT_WORLD));
PetscCall(VecView(left,PETSC_VIEWER_STDOUT_WORLD));

PetscCall(MatView(A,PETSC_VIEWER_DRAW_WORLD));

/*
//create
PetscInt length = 6;
PetscScalar val = 1.0/2.0;
Vec x;
PetscCall(VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,length,&x));
//PetscCall(VecView(x,PETSC_VIEWER_STDOUT_WORLD));

//get the range
PetscInt istart,iend;
PetscCall(VecGetOwnershipRange(x, &istart, &iend));
PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD,"s %d e %d  \n", istart, iend ));
//PetscCall(VecSet(x,val));

for (int i=istart; i<iend; i++)
PetscCall(VecSetValue(x,i,i,INSERT_VALUES));

PetscCall(VecAssemblyBegin(x));
PetscCall(VecAssemblyEnd(x));

//PetscCall(VecView(x,PETSC_VIEWER_STDOUT_WORLD));




PetscCall(VecDuplicate(x,&y));
PetscCall(VecCopy(x,y));
*/
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
/*
PetscCall(VecView(y,PETSC_VIEWER_STDOUT_WORLD));

//PetscScalar val;
//PetscCall(VecDot(x,y,&val));
//printf("dot prod:%lf",val);


PetscScalar product,norm;

PetscCall(VecDot(x,y,&product));
PetscCall(VecNorm(y,NORM_2,&norm));
PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Hello: %lf\t%lf\n",product,norm*norm));

*/


//PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Rank %d out of %d says hello \n", rank, size ));
PetscCall(PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT));

//PetscCall(VecDestroy())
PetscCall(PetscFinalize());

  return 0;
}
