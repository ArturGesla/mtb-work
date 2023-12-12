ALL: mat2

-include ../arch/make_inc

SOURCES    = mat2.F90
OBJ        = $(SOURCES:.F90=.o)
EXE        = mat2.exe
CLEANFILES = ${OBJ} ${EXE}

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

mat2: ${OBJ}
	${FLINKER} -o ${EXE} ${OBJ} ${PETSC_LIB}
