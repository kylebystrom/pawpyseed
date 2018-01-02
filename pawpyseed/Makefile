#!/bin/sh

projector:
	rm -f $(PWD)/pawpy.so
	source ~/.bashrcmkl
	icc -shared -o pawpy.so projector.c reader.c utils.c -L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -ldl -I${MKLROOT}/include -lm -fopenmp -std=c99 -O3 -fPIC -Wall

debug:
	rm -f $(PWD)/pawpy.so
	source ~/.bashrc2
	icc -shared -o pawpy.so projector.c reader.c utils.c -mkl=sequential -lm -std=c99 -O0 -fPIC -g -Wall

readtest:
	icc -o rtest readtest.c reader.c utils.c -lm -std=c99 -O0 -fPIC -g

tester:
	rm -f $(PWD)/pawpy.so
	gcc -shared -o pawpy.so reader.c utils.c -lm -std=c99 -O3 -fPIC
