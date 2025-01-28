CC = gcc
CFLAGS = -Wall -pedantic -std=gnu99 -Wextra -g -O0
DEBUG = -g
LDFLAGS = -lm -ffast-math

.DEFAULT_GOAL := ljmodel

.PHONY: debug clean

debug: CFLAGS += $(DEBUG)
debug: ljmodel
 
lennard_jones.o: lennard_jones.c
		$(CC) $(CFLAGS) -c $^ $(LDFLAGS) -o $@                                                  
                                                                             
ljmodel: lj_model_optimised.c lennard_jones.o               
		$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@                                 

clean:
		rm -f ljmodel *.o