#!/bin/sh

projector:
	rm -f $(PWD)/pawpy.so
	source ~/.bashrc2
	icc -c projector.c reader.c utils.c -lm -fopenmp -std=c99 -O3 -fPIC
	icc -shared -o pawpy.so projector.c reader.c utils.c -lm -fopenmp -std=c99 -O3 -fPIC

