ALL: dmda

-include ../../arch/make_inc

SOURCES    = dmda.F90
OBJ        = $(SOURCES:.F90=.o)
EXE        = dmda.exe
CLEANFILES = ${OBJ} ${EXE}

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

dmda: ${OBJ}
	${FLINKER} -o ${EXE} ${OBJ} ${PETSC_LIB}
