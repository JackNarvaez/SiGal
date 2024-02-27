CC=mpirun
PC=python3
SRC=galaxy.py
Np=4

setup: $(SRC)
	@$(CC) -np $(Np) $(PC) $(SRC)
	
clean:
	@rm -f ../Data/*.txt;\
	rm -f ../galaxy
