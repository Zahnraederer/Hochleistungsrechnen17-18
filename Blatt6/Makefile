.PHONY: clean

CC = mpicc
CFLAGS = -std=gnu11 -O3 -Wall -Werror -ggdb
LDFLAGS = $(CFLAGS)

all: timempi timempi2  

timempi2: timempi2.o
	$(CC) $(LDFLAGS) $@.o -o $@

timempi: timempi.o
	$(CC) $(LDFLAGS) $@.o -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean: 
	rm -f *.o *~ timempi timempi2
