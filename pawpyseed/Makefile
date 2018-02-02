#!/bin/sh

projector:
	rm -f $(PWD)/pawpy.so
	icc -shared -o pawpy.so projector.c fft.c quadrature.c reader.c utils.c radial.c tests.c \
		-std=c11 -L${MKLROOT}/lib/intel64 -lmkl_rt -fopenmp \
		-lpthread -ldl -I${MKLROOT}/include -lm -O3 -fPIC -Wall

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
	icc -o memtest memtest.c projector.c fft.c quadrature.c reader.c tests.c utils.c radial.c -mkl=sequential -L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -ldl -I${MKLROOT}/include -lm -std=c99 -O0 -fPIC -Wall -g

