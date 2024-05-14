from numpy import loadtxt, array, fromfile, float64, array
import matplotlib.pyplot as plt
import matplotlib.animation as animation
plt.style.use('dark_background')

def read_parameters(prmts):
    """------------------------------------------------------------------------
    Parameters
    ---------------------------------------------------------------------------
    Prmts:
    0       # Nb: Number of bodies in the galaxy
    1       # M : Mass of Galaxy
    2       # r : Radius of Galaxy
    3       # dt: Time step
    4       # steps: Evolution steps
    5       # jump: Data storage interval
    ------------------------------------------------------------------------"""
    data = loadtxt("../input")
    return data[prmts]

def read_data(file_path, N, coords):
    data = fromfile(file_path, dtype=float64)
    x = data[coords[0]:3*N+coords[0]:3]
    y = data[coords[1]:3*N+coords[1]:3]
    return x, y

N, R, dt, steps, jump = read_parameters([0, 2, 3, 4, 5])
N = int(N)
jump = int(jump)
steps = int(steps)//jump

#plot
xi = 0
yi = 1

# Extract the initial positions of the particles
initial_x, initial_y = read_data("../Data/Ev_0", N, [xi,yi])

def update_scatter(num, scatter, message):
    message.set_text(f" t = {(num+1)*jump*dt: .2f} yr")
    x, y = read_data("../Data/Ev_"+str(num), N, [xi,yi])
    scatter._offsets = array([x, y]).transpose()
    return scatter, message

# Attaching 3D axis to the figure
fig = plt.figure()
ax = fig.add_subplot(111)

message = ax.text(0.0, 0.9, "", transform=ax.transAxes)

# Creating a scatter plot for N bodies
scatter = ax.scatter(initial_x, initial_y, c='white', marker='.', s=1)
boundaries = 7.0 * R

# Setting the axes properties
if boundaries:
    ax.set_xlim([-boundaries, boundaries])
    ax.set_ylim([-boundaries, boundaries])
    ax.set_xlabel('x [au]')
    ax.set_ylabel('y [au]')

# Creating the Animation object
scatter_ani = animation.FuncAnimation(fig, update_scatter, frames=steps, fargs=(scatter, message), interval=1, blit=True, repeat=False)

scatter_ani.save('../Data/Evolution2D.gif', writer='pillow', fps=10)
plt.show()