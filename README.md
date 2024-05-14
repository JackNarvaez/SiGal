SiGal
=====
>Simulations of Galaxies using N Body Simulations. 

## Compilation

Clone this repository on your machine. Then, navigate to it and run the commands:

`$ make setup`

`$ mpirun -np <NP> ./_setup`

`$ make galaxy`

`$ mpirun -np <NP> ./galaxy`

where `<NP>` is the number of processes. To visualize the evolution, one can use the command `$ make animation3D`.
