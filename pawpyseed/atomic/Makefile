#!/bin/sh

ICC_LIBS = -L${MKLROOT}/lib/intel64
ICC_INCS = -I${MKLROOT}/include
ICC_FLAGS = -std=c11 -DMKL_Complex16="double complex" -lmkl_rt -fopenmp -lpthread -ldl -lm -O3 -fPIC -Wall
ICC_SRC = hf.c hf_utils.c condon_shortley.c
ICC_OBJ = hf.o hf_utils.o condon_shortley.o

GCC_LIBS = 
GCC_INCS = 
GCC_FLAGS = -std=c11 -fopenmp -lpthread -ldl -lm -O3 -fPIC -Wall

hfc: $(ICC_SRC)
	icc -shared -c $^ $(ICC_FLAGS) $(ICC_INCS) $(ICC_LIBS)
	icc -shared -o pawpy.so $(ICC_OBJ) $(ICC_FLAGS) $(ICC_INCS) $(ICC_LIBS)
