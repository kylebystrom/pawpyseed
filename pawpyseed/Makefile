#!/bin/sh

projector:
	rm -f $(PWD)/pawpy.so
	source ~/.bashrc
	icc -shared -o pawpy.so projector.c fft.c quadrature.c reader.c utils.c radial.c tests.c -std=c11 -L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -ldl -I${MKLROOT}/include -lm -O3 -fPIC -Wall

debug:
	rm -f $(PWD)/pawpy.so
	source ~/.bashrc2
	icc -shared -o pawpy.so projector.c reader.c utils.c -mkl=sequential -lm -std=c99 -O0 -fPIC -g -Wall

readtest:
	icc -o rtest readtest.c reader.c utils.c -lm -std=c99 -O0 -fPIC -g

tester:
	rm -f $(PWD)/pawpy.so
	icc -shared -o pawpy.so reader.c utils.c tests.c -lm -std=c99 -O3 -fPIC

forval:
	rm -f $(PWD)/memtest
	source ~/.bashrc
	icc -o memtest memtest.c projector.c fft.c quadrature.c reader.c tests.c utils.c radial.c -mkl=sequential -L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -ldl -I${MKLROOT}/include -lm -std=c99 -O0 -fPIC -Wall -g

formem:
	rm -f $(PWD)/pawpy.so
	source ~/.bashrcmkl
	icc -o memtest projector.c reader.c utils.c memtest.c -L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -ldl -I${MKLROOT}/include -lm -std=c99 -O0 -fPIC -Wall -g


