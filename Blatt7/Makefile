.PHONY: clean

CC = mpicc
CFLAGS = -std=c11 -O3 -Wall -Wextra -Werror -ggdb
LDFLAGS = $(CFLAGS)

all: circle  

circle: circle.o
	$(CC) $(LDFLAGS) $@.o -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean: 
	rm -f *.o *~ circle
