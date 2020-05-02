#
# Makefile for ct-cbn
#

VERSION = 0.1.00

#CC = gcc-4.2 -Wall -O3 -march=core2 -fopenmp  #-pg
CC = gcc -Wall -ggdb -mcmodel=large # -O3 -p #-pg


all: clonex

clonex.o: clonex.c
	$(CC)  -I/usr/local/include -c clonex.c

clonex: clonex.o 
	$(CC) clonex.o -L/usr/local/lib -lgsl -lgslcblas -lm -o $@

clean:
	rm -f a.out core *.o clonex *~ \#*



