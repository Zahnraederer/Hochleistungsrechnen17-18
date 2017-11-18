.PHONY: clean clean-all

# Common definitions
CC = gcc

# Compiler flags, paths and libraries
CFLAGS_SEQ = -std=c99 -pedantic -Wall -Wextra -O2 -ggdb
CFLAGS_OMP = -std=c99 -pedantic -Wall -Wextra -Ofast -flto -march=native -fopenmp
LFLAGS_SEQ = $(CFLAGS_SEQ)
LFLAGS_OMP = $(CFLAGS_OMP)
LIBS   = -lm

OBJS_SEQ = partdiff-seq.o askparams-seq.o
OBJS_OMP = partdiff-openmp.o askparams-omp.o

# Targets ...
all: partdiff-seq\
	partdiff-openmp-zeilen partdiff-openmp-spalten partdiff-openmp-element partdiff-openmp\
	partdiff-openmp 

partdiff-seq: partdiff-seq.o askparams-seq.o
	$(CC) $(LFLAGS_SEQ) -o $@ $^ $(LIBS)
	
partdiff-openmp: partdiff-openmp-zeilen.o askparams-omp.o
	$(CC) $(LFLAGS_OMP) -o $@ $^ $(LIBS)
	
partdiff-openmp-zeilen: partdiff-openmp-zeilen.o askparams-omp.o
	$(CC) $(LFLAGS_OMP) -o $@ $^ $(LIBS)

partdiff-openmp-spalten: partdiff-openmp-spalten.o askparams-omp.o
	$(CC) $(LFLAGS_OMP) -o $@ $^ $(LIBS)

partdiff-openmp-element: partdiff-openmp-element.o askparams-omp.o
	$(CC) $(LFLAGS_OMP) -o $@ $^ $(LIBS)

# Object files...
partdiff-seq.o: partdiff-seq.c
	$(CC) -c $(CFLAGS_SEQ) $<
	
partdiff-openmp-zeilen.o: partdiff-openmp-zeilen.c
	$(CC) -c $(CFLAGS_OMP) $<
	
partdiff-openmp-spalten.o: partdiff-openmp-spalten.c
	$(CC) -c $(CFLAGS_OMP) $<

partdiff-openmp-element.o: partdiff-openmp-element.c
	$(CC) -c $(CFLAGS_OMP) $<

askparams-seq.o: askparams.c
	$(CC) -c $(CFLAGS_SEQ) $< -o $@
	
askparams-omp.o: askparams.c
	$(CC) -c $(CFLAGS_OMP) $< -o $@

# PHONY...
clean:
	$(RM) *.o *~

clean-all: clean
	rm -f partdiff-seq partdiff-openmp-zeilen partdiff-openmp-spalten partdiff-openmp-element partdiff-openmp

