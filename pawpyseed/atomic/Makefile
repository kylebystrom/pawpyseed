#!/bin/sh

ICC_LIBS = -L${MKLROOT}/lib/intel64
ICC_INCS = -I${MKLROOT}/include
ICC_FLAGS = -std=c11 -lmkl_rt -fopenmp -lpthread -ldl -lm -O3 -fPIC -Wall
ICC_SRC = hf.c condon_shortley.c
ICC_OBJ = hf.o condon_shortley.o

GCC_LIBS = 
GCC_INCS = 
GCC_FLAGS = -std=c11 -fopenmp -lpthread -ldl -lm -O3 -fPIC -Wall

hfc:
	icc -shared -c $^ $(ICC_FLAGS) $(ICC_INCS) $(ICC_LIBS)
	icc -shared -o pawpy.so *.o $(ICC_FLAGS) $(ICC_INCS) $(ICC_LIBS)
