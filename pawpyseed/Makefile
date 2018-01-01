#!/bin/sh

projector:
	rm -f $(PWD)/pawpy.so
	source ~/.bashrc2
	icc -shared -o pawpy.so projector.c reader.c utils.c -lm -fopenmp -std=c99 -O3 -fPIC

debug:
	rm -f $(PWD)/pawpy.so
	source ~/.bashrc2
	icc -shared -o pawpy.so projector.c reader.c utils.c -lm -std=c99 -O0 -fPIC -g

tester:
	rm -f $(PWD)/pawpy.so
	gcc -shared -o pawpy.so reader.c utils.c -lm -std=c99 -O3 -fPIC
