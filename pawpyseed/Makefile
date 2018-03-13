#!/bin/sh

ICC_LIBS = -L${MKLROOT}/lib/intel64
ICC_INCS = -I${MKLROOT}/include
ICC_FLAGS = -std=c11 -lmkl_rt -fopenmp -lpthread -ldl -lm -O3 -fPIC -Wall
ICC_SRC = utils.c radial.c sbt.c reader.c quadrature.c fft.c pseudoprojector.c projector.c tests.c
ICC_OBJ = utils.o radial.o sbt.o reader.o quadrature.o fft.o pseudoprojector.o projector.o tests.o

GCC_LIBS = 
GCC_INCS = 

projector: $(ICC_SRC)
	rm -f $(PWD)/pawpy.so
	icc -shared -c $@ $^ $(ICC_FLAGS) $(ICC_INCS) $(ICC_LIBS)
	icc -shared -o pawpy.so *.o $(ICC_FLAGS) $(ICC_INCS) $(ICC_LIBS)

projector_gcc:
	rm -f $(PWD)/pawpy.so
	gcc -shared -o pawpy.so pseudoprojector.c projector.c gsl_fft.c quadrature.c \
		reader.c utils.c radial.c tests.c sbt.c \
		-std=c11 -fopenmp \
		-lpthread -ldl -lm -O3 -fPIC -Wall

lite:
	rm -f $(PWD)/pawpy.so
	icc -shared -o pawpy.so pseudoprojector.c utils.c -std=c11 -fopenmp -lm -O3 -Wall -fPIC

lite_gcc:
	rm -f $(PWD)/pawpy.so
	gcc -shared -o pawpy.so pseudoprojector.c utils.c -std=c11 -fopenmp -lm -O3 -Wall -fPIC

test:
	rm -f $(PWD)/pawpy.so
	icc -shared -o pawpy.so reader.c utils.c tests.c -lm -std=c99 -O3 -fPIC

memtest:
	rm -f $(PWD)/memtest
	icc -o memtest memtest.c pseudoprojector.c projector.c fft.c quadrature.c reader.c tests.c utils.c radial.c sbt.c -mkl=sequential -L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -ldl -I${MKLROOT}/include -lm -std=c99 -O0 -fPIC -Wall -g

