CC=mpic++
CFLAGS=-Wall -Werror
COPTS=-O3
DEPS=nbodies.h structures.h tree.h data.h parallel.h
SRC=main.c nbodies.c tree.c data.c parallel.c
OBJ=${SRC:.c=.o}

%.o: %.c $(DEPS)
	@$(CC) $(CFLAGS) $(COPTS) -o $@ -c $<

galaxy: $(OBJ)
	@$(CC) $(CFLAGS) $(COPTS) -o ../$@ $(OBJ)

clean:
	@rm -f *.o *.mod *.out