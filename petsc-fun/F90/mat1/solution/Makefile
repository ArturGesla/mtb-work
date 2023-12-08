ALL: mat1

-include ../../arch/make_inc

SOURCES    = mat1.F90
OBJ        = $(SOURCES:.F90=.o)
EXE        = mat1.exe
CLEANFILES = ${OBJ} ${EXE}

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

mat1: ${OBJ}
	${FLINKER} -o ${EXE} ${OBJ} ${PETSC_LIB}
