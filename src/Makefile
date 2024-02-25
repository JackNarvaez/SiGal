CC=mpic++
CFLAGS=-Wall -Werror
COPTS=-O0
DEPS=NBodies.h
SRC=main.c NBodies.c
OBJ=${SRC:.c=.o}

%.o: %.c $(DEPS)
	@$(CC) $(CFLAGS) $(COPTS) -o $@ -c $<

galaxy: $(OBJ)
	@$(CC) $(CFLAGS) $(COPTS) -o ../$@ $(OBJ)

clean:
	@rm -f *.o *.mod *.out