Nb=1000		# Number of bodies in the galaxy
Np=4           # Number of processes
i=0.1 		    # Inclination of the galaxy (rads/pi)
w=0.6 		    # Angle in the xy plane of the galaxy (rads/pi)
steps=10000 	# Evolution steps
jump=100		# Data storage interval
dt=0.001		# Time step
rad=5000		# Radius of Galaxy (AU)

cd setup
mpirun -n ${Np} python3 galaxy.py ${Nb} ${i} ${w} ${rad}
cd ../src
mpic++ main.c NBodies.c -o Galaxy.x
mpirun -np ${Np} ./Galaxy.x ${steps} ${dt} ${jump} ${Nb}
cd ../plot
conda activate
python3 Animation.py ${Nb} ${dt} ${jump} ${rad}