#run make GSL=0 for compilation without GSL
ALL: build

SRC = src/matrix_constructions.c \
			src/create_index.c \
			src/main.c \
			src/generate_parameters.c \

-include ${PETSC_DIR}/lib/petsc/conf/variables

MAINFILE		:= hdiv
CPPFLAGS    := -I./include ${PETSC_CCPPFLAGS}
CFLAGS	    := ${CC_FLAGS}
OBJ         := $(SRC:.c=.o)
LIB 	      := ${PETSC_TAO_LIB}


ifneq ($(GSL),0)
	CFLAGS += -DHAVE_GSL
	LIB += -lgsl -lgslcblas
endif

src/%.o: %.c
	@echo ${PETSC_COMPILE_SINGLE}

${MAINFILE}: ${OBJ}
	-${CLINKER} -o ${MAINFILE} ${OBJ} ${LIB} 

build: ${MAINFILE}

clean:
	-@${RM} ${OBJ}
	-@${RM} ${OBJ:.o=.d}
	-@${RM} ${MAINFILE}

.PHONY: ALL build clean

