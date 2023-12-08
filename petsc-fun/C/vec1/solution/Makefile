ALL: vec1

-include ../../arch/make_inc

SOURCES    = vec1.c
OBJ        = $(SOURCES:.c=.o)
EXE        = vec1.exe
CLEANFILES = ${OBJ} ${EXE}

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

vec1: ${OBJ}
	${CLINKER} -o ${EXE} ${OBJ} ${PETSC_LIB}
