CC=mpic++
CFLAGS=-Wall -Werror
COPTS=-O0
DEPS=hernquist.h ../src/structures.h ../src/data.h
SRC=main.c hernquist.c ../src/data.c
OBJ=${SRC:.c=.o}

%.o: %.c $(DEPS)
	@$(CC) $(CFLAGS) $(COPTS) -o $@ -c $<

setup: $(OBJ)
	@$(CC) $(CFLAGS) $(COPTS) -o ../_$@ $(OBJ)

clean:
	@rm -f ../Data/*.txt *.o;\
	rm -f ../galaxy ../_setup