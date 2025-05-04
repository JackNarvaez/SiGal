CC=mpicc
CFLAGS=-Wall -Werror
CFLAGSPROF=-g -pg -Wall -Werror
COPTS=-O3
LDFLAGS=-L/usr/local/lib -lgsl -lgslcblas -lm
DEPS=nbodies.h structures.h tree.h data.h parallel.h
SRC=main.c nbodies.c tree.c data.c parallel.c
OBJ=${SRC:.c=.o}

%.o: %.c $(DEPS)
	@$(CC) $(CFLAGS) $(COPTS) -o $@ -c $<

galaxy: $(OBJ)
	@$(CC) $(CFLAGS) $(COPTS) -o ../$@ $(OBJ) $(LDFLAGS)

data.o: data.c $(DEPS)
	$(CC) $(CFLAGS) $(COPTS) -o $@ -c $

clean:
	@rm -f *.o *.mod *.out