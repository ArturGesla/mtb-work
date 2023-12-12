ALL: solver

-include ../arch/make_inc

SOURCES    = solver.F90
OBJ        = $(SOURCES:.F90=.o)
EXE        = solver.exe
CLEANFILES = ${OBJ} ${EXE}

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

solver: ${OBJ}
	${FLINKER} -o ${EXE} ${OBJ} ${PETSC_LIB}
