ALL: mat2

-include ../arch/make_inc

SOURCES    = mat2.c
OBJ        = $(SOURCES:.c=.o)
EXE        = mat2.exe
CLEANFILES = ${OBJ} ${EXE}

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

mat2: ${OBJ}
	${CLINKER} -o ${EXE} ${OBJ} ${PETSC_LIB}
