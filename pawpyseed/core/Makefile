#!/bin/sh

LIBS = -L${MKLROOT}/lib/intel64
INCS = -I${MKLROOT}/include
FLAGS = -std=c11 -lmkl_rt -fopenmp -lpthread -ldl -lm -O3 -fPIC -Wall -DMKL_Complex16="double complex" -DMKL_Complex8="float complex"
SRC = utils.c gaunt.c radial.c sbt.c reader.c quadrature.c linalg.c density.c pseudoprojector.c projector.c
OBJ = utils.o gaunt.o radial.o sbt.o reader.o quadrature.o linalg.o density.c pseudoprojector.o projector.o
TST_FLAGS = -std=c11 -lmkl_rt -lpthread -ldl -lm -O3 -fPIC -Wall -DMKL_Complex16="double complex" -DMKL_Complex8="float complex" 
TST_SRC = utils.c gaunt.c radial.c sbt.c reader.c quadrature.c linalg.c density.c pseudoprojector.c projector.c tests.c
TST_OBJ = utils.o gaunt.o radial.o sbt.o reader.o quadrature.o linalg.o density.c pseudoprojector.o projector.o tests.o

pawpyinst:
	$(PAWPYCC) -shared -c $(SRC) $(FLAGS) $(INCS) $(LIBS)
	$(PAWPYCC) -shared -o pawpy.so $(OBJ) $(FLAGS) $(INCS) $(LIBS)

testsuite:
	$(PAWPYCC) -shared -c $(TST_SRC) $(TST_FLAGS) $(INCS) $(LIBS)
	$(PAWPYCC) -shared -o pawpy.so $(TST_OBJ) $(TST_FLAGS) $(INCS) $(LIBS)

mem:
	$(PAWPYCC) -o memtest memtest.c gaunt.c pseudoprojector.c projector.c \
	linalg.c quadrature.c reader.c tests.c utils.c radial.c sbt.c \
	-L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread \
	-ldl -I${MKLROOT}/include -lm -std=c99 -O0 -fPIC -Wall -g \
	-DMKL_Complex16="double complex" -DMKL_Complex8="float complex"

lite:
	rm -f $(PWD)/pawpy.so
	$(PAWPYCC) -shared -o lpawpy.so pseudoprojector.c utils.c -std=c11 -fopenmp -lm -O3 -Wall -fPIC

#test:
#	rm -f $(PWD)/pawpy.so
#	$(PAWPYCC) -shared -o otpawpy.so reader.c utils.c tests.c -lm -std=c99 -O3 -fPIC
