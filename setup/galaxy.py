from numpy import array, zeros, random, sqrt, exp, pi, sin, cos, arctan, loadtxt
import scipy.integrate as integrate
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from mpi4py import MPI
from sys import argv

def func(x,Distribution,Point): 
    """------------------------------------------------------------------------
    func: 
    Equation that follows the point of the wanted distribution that matches the 
    random one of a uniform distribution
    ---------------------------------------------------------------------------
    Arguments:
       x            : Random variable in the wanted distribution (unkonwn)
       Distribution : Wanted distribution
       Point        : Random variable in the uniform distribution
    ------------------------------------------------------------------------"""
    return integrate.quad(Distribution,0,x)[0]-Point

def spiral_galaxy(N, max_mass, BHM, center, ini_radius, beta, alpha):
    '''-----------------------------------------------------------------------
    spiral_galaxy:
    Use a radial distrubution of masses proportional to the brightness surface
    distributation to create a plain Bulb and Disk resembling an spiral galaxy
    --------------------------------------------------------------------------
    Arguments:
       N            : Number of particles
       max_mass     : Biggest mass of the stars in the system 
       BHM          : Black Hole's mass
       center       : Black Hole position
       ini_radius   : Galaxy radius
       beta         : Inclination
       alpha        : Angle in the x,y plane
    ------------------------------------------------------------------------'''
    # Generates N random particles 
    positions = zeros(3*N)
    velocity= zeros(3*N)
    # Random masses varies between 1 solar mass and max_mass solar masses
    masses = random.random(N)*(max_mass-1.) + 1.
    #Parameters of the model of density of starts
    initial_density=.1
    const_bulb=.3
    const_disc=.8
    bulb_radius=0.2
    #Model of density normalized
    f1 = lambda x: initial_density*exp(-x**(1/4)/const_bulb)        #Bulge
    f2 = lambda x: f1(bulb_radius)*exp(-(x-bulb_radius)/const_disc) #Disc
    f = lambda x:  f1(x) if x<bulb_radius else f2(x)                #Piecewise 
    norm = integrate.quad(f,0,1)[0]                                  
    uf=lambda x: f(x)/norm                                          #Density function with integral=1
    #Random angle generation
    gamma = random.random(N)*2*pi
    #Random width
    width = .05*ini_radius                                          #Half of with in relation to the radius of the galaxy
    gross  = random.random(N)*2*width-width
    temp = beta
    #Uniform distribution to get random points
    Uniform = random.random(N)
    #Empty array for the points mapped from the uniform distribution
    Map=zeros(N)   
    for i in range(N):
        #Calls the function that maps the ramdom points to the wanted distribution for the radius 
        Map[i]=fsolve(func,0,args=(uf,Uniform[i]))*ini_radius
        #In case radius is to small, add a value to avoid particles scaping
        if Map[i] < 50:
            Map[i] += 50
        #Creates an elipsoid in the region of the bulge
        if Map[i] < bulb_radius*ini_radius:
            a = 0.18*ini_radius
            bulg_countour = a*sqrt(1-(Map[i]/(bulb_radius*ini_radius))**2)
            gross[i] = random.random(1)*2*bulg_countour-bulg_countour
        #Adjustment for width
        beta += arctan(gross[i]/Map[i])
        Map[i] = sqrt(Map[i]**2+gross[i]**2)
        #Change to cartesian coordinates
        positions[3*i+0] = Map[i]*(cos(gamma[i])*cos(alpha)+
                                  sin(gamma[i])*cos(beta)*sin(alpha)) + center[0]
        positions[3*i+1] = Map[i]*(sin(gamma[i])*cos(beta)*cos(alpha)-
                                  cos(gamma[i])*sin(alpha))+ center[1]
        positions[3*i+2] = Map[i]*sin(gamma[i])*sin(beta) + center[2]
        # Keplerina velocity in the plain of the disc 
        Kep_v = sqrt(G*BHM/Map[i])
        vec_vel=array([-Map[i]*(sin(gamma[i])*cos(alpha)-cos(gamma[i])*cos(beta)*sin(alpha)),
                       Map[i]*(cos(gamma[i])*cos(beta)*cos(alpha)+sin(gamma[i])*sin(alpha)), 
                       Map[i]*cos(gamma[i])*sin(beta)])/Map[i]
        velocity[3*i+0] = Kep_v*vec_vel[0]
        velocity[3*i+1] = Kep_v*vec_vel[1]
        velocity[3*i+2] = Kep_v*vec_vel[2]
        beta = temp
    return masses, positions, velocity

def Data(pId, size, max_mass, BHM, BHposition, ini_radius, beta, alpha):
    """------------------------------------------------------------------------
    Data: 
    Creates the information of each particle and saves into a file
    ---------------------------------------------------------------------------
    Arguments:
       pId          : Process id
       size         : Number of particles to each process
       max_mass     : Biggest mass of the stars in the system 
       BHM          : Black Hole's mass
       BHposition   : Black Hole position
       ini_radius   : Galaxy radius
       beta         : Inclination
       alpha        : Angle in the x,y plane
    ------------------------------------------------------------------------"""
    data = open("./../Data/data"+str(pId)+".txt", "w")
    data.write("#\t"+str(size)+"\n")
    if pId == 0:
        data.write(str(BHposition[0])+"\t"+str(BHposition[1])+"\t"+str(BHposition[2])+"\t0\t0\t0\t"+str(BHM)+"\n")
        size -= 1
    masses, positions, velocity = spiral_galaxy(size, max_mass, BHM, BHposition, ini_radius, beta, alpha)
    #Save
    for ii in range(size):
        data.write(str("{:.6f}".format(positions[3*ii+0]))+"\t"+
                   str("{:.6f}".format(positions[3*ii+1]))+"\t"+
                   str("{:.6f}".format(positions[3*ii+2]))+"\t"+
                   str("{:.6f}".format(velocity[3*ii+0]))+"\t"+
                   str("{:.6f}".format(velocity[3*ii+1]))+"\t"+
                   str("{:.6f}".format(velocity[3*ii+2]))+"\t"+
                   str("{:.6f}".format(masses[ii]))+"\n")
    data.close()

def read_parameters(prmts):
    """------------------------------------------------------------------------
    Parameters
    ---------------------------------------------------------------------------
    Prmts:
    0       # Nb: Number of bodies in the galaxy
    1       # i : Inclination of the galaxy (rads/pi)
    2       # w : Angle in the xy plane of the galaxy (rads/pi)
    3       # dt: Time step
    4       # r : Radius of Galaxy (AU)
    5       # steps: Evolution steps
    6       # jump: Data storage interval
    ------------------------------------------------------------------------"""
    data = loadtxt("../input")
    return data[prmts]


# ----------------------------MAIN-------------------------------- #
#MPI Initialize 
comm = MPI.COMM_WORLD
pId = comm.Get_rank()
nP = comm.Get_size()
# Gravitational constant in units of au^3 M_sun^-1 yr^-2
G = 4*pi**2
# Mass of the N bodies.
max_mass = 50. # Solar masses
# Parameters of the galaxy plane orientation 
N, i, w, ini_radius = read_parameters([0, 1, 2, 4])
beta=pi*i     #Inclination
alpha=pi*w    #Angle in the plain x,y
# Supermassive Central Black Hole data
BHM = 1.e6 # Solar masses
BHposition = array([0., 0., 0.]) # Location of the SBH
# Random Seed
random.seed(pId)
# Partition between processes
teil = float(N)/nP
size = int(teil*(pId+1))-int(teil*pId)
Data(pId, size, max_mass, BHM, BHposition, ini_radius, beta, alpha)
