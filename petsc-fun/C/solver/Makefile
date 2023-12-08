ALL: solver

-include ../arch/make_inc

SOURCES    = solver.c
OBJ        = $(SOURCES:.c=.o)
EXE        = solver.exe
CLEANFILES = ${OBJ} ${EXE}

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

solver: ${OBJ}
	${CLINKER} -o ${EXE} ${OBJ} ${PETSC_LIB}
