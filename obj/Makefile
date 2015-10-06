VPATH=../src/

CC=gcc#${HOME}/Programs/openmpi-1.8.3/build/local/bin/mpicc
CFLAGS=-O4 -g -Wall
CPPFLAGS=#-DMPI
MPIINCLUDE=#-I${HOME}/Programs/openmpi-1.8.3/build/local/include
INCLUDE=$(MPIINCLUDE)
MPILIBS=#-L${HOME}/Programs/openmpi-1.8.3/build/local/lib -lmpi
LAPACKLIBS=-llapack -lblas#-L${HOME}/Programs/lapack-3.5.0 -llapack -lblas
FFTWLIBS=-L${HOME}/Programs/fftw-3.3.4/build/local/lib -lfftw3
LIBS=$(MPILIBS) $(LAPACKLIBS) $(FFTWLIBS) -lc -lm -lgfortran

OBJS=sys.o constants.o numeric.o matrix.o randnum.o vector.o fourier.o lapack.o \
     parameters.o conjgrad.o davidson.o dmrg.o tddmrg.o
     
MAIN=fermichain.o

all: fermichain
clean:
	rm -f $(OBJS) $(MAIN) fermichainexec

fermichain: $(OBJS) $(MAIN)
	$(CC) $(CFLAGS) $(CPPFLAGS) $(INCLUDE) $(OBJS) $(MAIN) -o fermichainexec $(LIBS)
