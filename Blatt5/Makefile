# Common definitions
CC = gcc

# Compiler flags, paths and libraries
CFLAGS = -std=c11 -pedantic -Wall -Wextra -O3 -ggdb
CFLAGS_PTHREAD = -std=c11  -pedantic -Wall -Wextra -ggdb -lpthread
LFLAGS = $(CFLAGS)
LFLAGS_PTHREAD = $(CFLAGS_PTHREAD)
LIBS   = -lm

TGTS = partdiff-seq
OBJS = partdiff-seq.o askparams.o
TBJS = partdiff-pthread.o askparams.o


# Targets ...
all: partdiff-seq\
	partdiff-pthread

partdiff-seq: $(OBJS) Makefile
	$(CC) $(LFLAGS) -o $@ $(OBJS) $(LIBS)

partdiff-pthread: $(TBJS) Makefile
	$(CC) $(LFLAGS_PTHREAD) -o $@ $(TBJS) $(LIBS)

partdiff-pthread.o: partdiff-pthread.c Makefile

partdiff-seq.o: partdiff-seq.c Makefile

askparams.o: askparams.c Makefile

# Rule to create *.o from *.c
%.o: %.c
	$(CC) -c $(CFLAGS) $*.c

clean:
	$(RM) $(OBJS)
	$(RM) $(TGTS)
