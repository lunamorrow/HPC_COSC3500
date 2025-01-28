MPI_INC = /usr/include/openmpi-x86_64
MPI_LIB = /usr/lib64/openmpi/lib

CC = mpicc
CFLAGS = -Wall -std=gnu99 -g -O3 -I$(MPI_INC)
DEBUG = -g
LDFLAGS = -lm -ffast-math -fopenmp -L$(MPI_LIB) 

.DEFAULT_GOAL := ljmodel

.PHONY: debug clean

debug: CFLAGS += $(DEBUG)
debug: ljmodel
 
lennard_jones.o: lennard_jones.c
		$(CC) $(CFLAGS) -c $^ $(LDFLAGS) -o $@                                                  									 
													                         
ljmodel: lj_model.c lennard_jones.o                 
		$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@ 

clean:
		rm -f ljmodel *.o