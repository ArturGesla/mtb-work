ALL: dmda

-include ../arch/make_inc

SOURCES    = dmda.c
OBJ        = $(SOURCES:.c=.o)
EXE        = dmda.exe
CLEANFILES = ${OBJ} ${EXE}

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

dmda: ${OBJ}
	${CLINKER} -o ${EXE} ${OBJ} ${PETSC_LIB}
