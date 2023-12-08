ALL: vec2a vec2b vec2c

-include ../../arch/make_inc

SOURCESA   = vec2a.F90
OBJA       = $(SOURCESA:.F90=.o)
EXEA       = vec2a.exe
SOURCESB   = vec2b.F90
OBJB       = $(SOURCESB:.F90=.o)
EXEB       = vec2b.exe
SOURCESC   = vec2c.F90
OBJC       = $(SOURCESC:.F90=.o)
EXEC       = vec2c.exe
CLEANFILES = ${OBJA} ${EXEA} ${OBJB} ${EXEB} ${OBJC} ${EXEC}

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

vec2a: ${OBJA}
	${FLINKER} -o ${EXEA} ${OBJA} ${PETSC_LIB}

vec2b: ${OBJB}
	${FLINKER} -o ${EXEB} ${OBJB} ${PETSC_LIB}

vec2c: ${OBJC}
	${FLINKER} -o ${EXEC} ${OBJC} ${PETSC_LIB}
