CC=mpic++
CFLAGS=-Wall -Werror
COPTS=-O0
DEPS = ../src/structures.h ../src/data.h utils.h plummer.h hernquist.h exponential.h kepler.h Kuzmin.h merger.h
SRC=main.c ../src/data.c utils.c plummer.c hernquist.c exponential.c kepler.c Kuzmin.c merger.c
OBJ=${SRC:.c=.o}

%.o: %.c $(DEPS)
	@$(CC) $(CFLAGS) $(COPTS) -o $@ -c $<

setup: $(OBJ)
	@$(CC) $(CFLAGS) $(COPTS) -o ../_$@ $(OBJ)

clean:
	@rm -f ../Data/*.txt *.o;\
	rm -f ../galaxy ../_setup