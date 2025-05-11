ALL: vec1

-include ../../arch/make_inc

SOURCES    = vec1.f90
OBJ        = $(SOURCES:.f90=.o)
EXE        = vec1.exe
CLEANFILES = ${OBJ} ${EXE}

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

vec1: ${OBJ}
	${FLINKER} -o ${EXE} ${OBJ} ${PETSC_LIB}
