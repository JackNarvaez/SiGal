CC=mpic++
CFLAGS=-Wall -Werror
LDFLAGS=-L/usr/local/lib -lgsl -lgslcblas -lm
COPTS=-O0
DEPS = ../src/structures.h ../src/data.h utils.h plummer.h hernquist.h exponential.h kepler.h kuzmin.h miyamoto.h merger.h mwg1.h
SRC=main.c ../src/data.c utils.c plummer.c hernquist.c exponential.c kepler.c kuzmin.c miyamoto.c merger.c mwg1.c
OBJ=${SRC:.c=.o}

%.o: %.c $(DEPS)
	@$(CC) $(CFLAGS) $(COPTS) -o $@ -c $<

setup: $(OBJ)
	@$(CC) $(CFLAGS) $(COPTS) -o ../_$@ $(OBJ) $(LDFLAGS)

clean:
	@rm -f ../Data/*.txt *.o;\
	rm -f ../galaxy ../_setup