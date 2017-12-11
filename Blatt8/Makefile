.PHONY: clean

CC = mpicc
CFLAGS = -std=c11 -O3 -Wall -Wextra -Werror -ggdb
LDFLAGS = $(CFLAGS)

all: partdiff-par

partdiff-par: partdiff-par.o displaymatrix-mpi.o
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

partdiff-par.o: partdiff-par.c
	$(CC) -c $(CFLAGS) $<

displaymatrix-mpi.o: displaymatrix-mpi.c
	$(CC) -c $(CFLAGS) $<

clean: 
	rm -f *.o *~ partdiff-par
