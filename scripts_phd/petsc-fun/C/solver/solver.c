#include <petsc.h>
int main(int argc, char **argv)
{



PetscMPIInt rank, size;
PetscFunctionBeginUser;
PetscCall(PetscInitialize(&argc,&argv,NULL,NULL));
PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));

//logging options
PetscLogStage stage1, stage2, stage3;
PetscLogStageRegister("Fill and assembly", &stage1);
PetscLogStageRegister("Solve", &stage2);
PetscLogStageRegister("Norm calc", &stage3);




PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Code to solve the linear system with distributed memory with PETSc.\n"));

PetscInt size1d;
PetscCall(PetscOptionsGetInt(NULL,NULL, "-size",&size1d, NULL));
PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD,"size1d:%d.\n",size1d));



//matrix definition - this would have to be changed
Mat A;

PetscInt n1d=size1d;
PetscInt n=n1d*n1d;
PetscScalar l=1;
PetscScalar h=l/n1d;

PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Parms: \t n: %d \t l: %lf \t h: %lf \n",n,l,h ));

PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
PetscCall(MatSetType(A,MATMPIAIJ));
PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n));

PetscInt istart,iend;
PetscCall(MatGetOwnershipRange(A,&istart,&iend));
PetscCall(PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Rank: %d \t size: %d \t istart: %d \t iend: %d\n",rank,size,istart, iend ));

//...
PetscLogStagePush(stage1);


PetscInt nloc=iend-istart;
PetscScalar center=4.0/h/h;;
PetscScalar offc=-1.0/h/h;

//memory is divided by lines
for (int i=istart; i<iend; i++)
{
//xi[i]=istart+i;
//vals[i]=1;
//printf("rank:%d\t el: %i\n",rank,i);
PetscCall(MatSetValues(A,1,&i,1,&i,&center,INSERT_VALUES));

if(i<n-1 && (i%n1d)!=(n1d-1))
{
int ip1=(i+1);
PetscCall(MatSetValues(A,1,&i,1,&ip1,&offc,INSERT_VALUES));
}
if(i<n-n1d)
{
int ip1=(i+n1d);
PetscCall(MatSetValues(A,1,&i,1,&ip1,&offc,INSERT_VALUES));
}
if(i>0 && (i%n1d)!=0)
{
int ip1=(i-1);
PetscCall(MatSetValues(A,1,&i,1,&ip1,&offc,INSERT_VALUES));
}
if(i>n1d-1)
{
int ip1=(i-n1d);
PetscCall(MatSetValues(A,1,&i,1,&ip1,&offc,INSERT_VALUES));
}

}


PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));

PetscCall(MatSetUp(A));

//PetscCall(MatView(A,PETSC_VIEWER_STDOUT_WORLD));
//PetscCall(MatView(A,PETSC_VIEWER_DRAW_WORLD));


//solver
KSP ksp;
KSPCreate(PETSC_COMM_WORLD,&ksp);
KSPSetOperators(ksp,A,A);
PetscCall(KSPSetFromOptions(ksp));
PetscCall(KSPSetUp(ksp));

Vec b,x;

PetscCall(VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,n,&b));
PetscCall(VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,n,&x));
//PetscCall(VecView(x,PETSC_VIEWER_STDOUT_WORLD));

for (int i=istart; i<iend; i++)
PetscCall(VecSetValue(b,i,1.0,INSERT_VALUES));

PetscCall(VecAssemblyBegin(b));
PetscCall(VecAssemblyEnd(b));

//PetscCall(VecView(b,PETSC_VIEWER_STDOUT_WORLD));

PetscLogStagePop();


PetscCall(KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD));

PetscLogStagePush(stage2);
PetscCall(KSPSolve(ksp,b,x));
PetscLogStagePop();



//PetscCall(VecView(x,PETSC_VIEWER_STDOUT_WORLD));
//PetscCall(VecView(x,PETSC_VIEWER_STDOUT_WORLD));


Vec Axmb;
PetscCall(VecDuplicate(x,&Axmb));
PetscCall(VecScale(b,-1));
PetscCall(MatMultAdd(A,x,b,Axmb));


PetscLogStagePush(stage3);

//norm
PetscScalar norm;
PetscCall(VecNorm(Axmb,NORM_2,&norm));
PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Hello: norm KSP Ax-b:%4.2e\n",norm));

PetscLogStagePop();






/*
//post proc
Mat X;
PetscCall(MatCreate(PETSC_COMM_WORLD,&X));

PetscCall(MatSetType(X,MATMPIDENSE));
PetscCall(MatSetSizes(X,PETSC_DECIDE,PETSC_DECIDE,n1d,n1d));

for (int i=0; i<n1d; i++){
for (int j=0; j<n1d; j++){
PetscScalar val;
PetscInt ind=i+j*n1d;

PetscCall(VecGetValues(x,1,&ind,&val));
PetscCall(MatSetValues(X,1,&i,1,&j,&val,INSERT_VALUES));
}}



PetscCall(MatAssemblyBegin(X,MAT_FINAL_ASSEMBLY));
PetscCall(MatAssemblyEnd(X,MAT_FINAL_ASSEMBLY));

PetscCall(MatSetUp(X));
PetscCall(MatView(X,PETSC_VIEWER_STDOUT_WORLD));
*/




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
