ALL: vec2a vec2b vec2c

-include ../arch/make_inc

SOURCESA   = vec2a.c
OBJA       = $(SOURCESA:.c=.o)
EXEA       = vec2a.exe
SOURCESB   = vec2b.c
OBJB       = $(SOURCESB:.c=.o)
EXEB       = vec2b.exe
SOURCESC   = vec2c.c
OBJC       = $(SOURCESC:.c=.o)
EXEC       = vec2c.exe
CLEANFILES = ${OBJA} ${EXEA} ${OBJB} ${EXEB} ${OBJC} ${EXEC}

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

vec2a: ${OBJA}
	${CLINKER} -o ${EXEA} ${OBJA} ${PETSC_LIB}

vec2b: ${OBJB}
	${CLINKER} -o ${EXEB} ${OBJB} ${PETSC_LIB}

vec2c: ${OBJC}
	${CLINKER} -o ${EXEC} ${OBJC} ${PETSC_LIB}
