.PHONY: clean

CC = mpicc
CFLAGS = -std=c11 -O3 -Wall -Wextra -Werror -ggdb -lm
LDFLAGS = $(CFLAGS)

all: partdiff-par

partdiff-par: partdiff-par.o askparams.o partdiff-seq.o #displaymatrix-mpi.o
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

partdiff-par.o: partdiff-par.c
	$(CC) -c $(CFLAGS) $<

displaymatrix-mpi.o: displaymatrix-mpi.c
	$(CC) -c $(CFLAGS) $<

askparams.o: askparams.c
	$(CC) -c $(CFLAGS) $<

partdiff-seq.o: partdiff-seq.c
	$(CC) -c $(CFLAGS) $<

clean: 
	rm -f *.o *~ partdiff-par
