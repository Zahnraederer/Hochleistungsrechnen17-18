# Common definitions
CC = gcc

# Compiler flags, paths and libraries
CFLAGS_SEQ = -std=c99 -pedantic -Wall -Wextra -O2 -ggdb
CFLAGS_POSIX = -std=c99 -pedantic -Wall -Wextra -Ofast -flto -march=native -lpthread
LFLAGS_SEQ = $(CFLAGS_SEQ)
LFLAGS_POSIX = $(CFLAGS_POSIX)
LIBS   = -lm

OBJS_SEQ = partdiff-seq.o askparams-seq.o
OBJS_POSIX = partdiff-posix.o askparams-posix.o

# Targets ...
all: partdiff-seq partdiff-posix

partdiff-seq: partdiff-seq.o askparams-seq.o
	$(CC) $(LFLAGS_SEQ) -o $@ $^ $(LIBS)
	
partdiff-posix: partdiff-posix.o askparams-posix.o
	$(CC) $(LFLAGS_POSIX) -o $@ $^ $(LIBS)

# Object files...
partdiff-seq.o: partdiff-seq.c
	$(CC) -c $(CFLAGS_SEQ) $<
	
partdiff-posix.o: partdiff-posix.c
	$(CC) -c $(CFLAGS_POSIX) $<
	
askparams-seq.o: askparams.c
	$(CC) -c $(CFLAGS_SEQ) $< -o $@
	
askparams-posix.o: askparams.c
	$(CC) -c $(CFLAGS_POSIX) $< -o $@

# PHONY...
clean:
	$(RM) *.o *~

clean-all: clean
	rm -f partdiff-seq partdiff-posix

