ALL: mat1

-include ../arch/make_inc

SOURCES    = mat1.c
OBJ        = $(SOURCES:.c=.o)
EXE        = mat1.exe
CLEANFILES = ${OBJ} ${EXE}

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

mat1: ${OBJ}
	${CLINKER} -o ${EXE} ${OBJ} ${PETSC_LIB}
