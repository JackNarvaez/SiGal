from numpy import loadtxt, array, transpose
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


N, R, dt, jump= read_parameters([0, 2, 3, 5])
N = int(N)
jump = int(jump)

x, y, m = loadtxt("../Data/Evolution.txt", usecols=[0, 1, 3], unpack=1)    # Evolution
steps = x.size // N # Printed steps

# Extract the initial positions of the particles
initial_x = x[:N]
initial_y = y[:N]

def update_scatter(num, scatter, message):
    message.set_text(f" t = {num*jump*dt: .2f} yr")
    scatter._offsets = (array([x[num*N:(num+1)*N], y[num*N:(num+1)*N]]).transpose())
    return scatter, message


# Attaching 3D axis to the figure
fig = plt.figure()
ax = fig.add_subplot(111)

message = ax.text(0.0, 0.9, "", transform=ax.transAxes)

# Creating a scatter plot for N bodies
scatter = ax.scatter(initial_x, initial_y, c='white', marker='.', s=1)
boundaries = 7.0 * R
ax.set_xlim([-boundaries, boundaries])
ax.set_ylim([-boundaries, boundaries])
plt.savefig("parallel.png")

# Setting the axes properties
if boundaries:
    ax.set_xlim([-boundaries, boundaries])
    ax.set_ylim([-boundaries, boundaries])
    ax.set_xlabel('x [au]')
    ax.set_ylabel('y [au]')

# Creating the Animation object
scatter_ani = animation.FuncAnimation(fig, update_scatter, frames=steps, fargs=(scatter, message),
                                      interval=1, blit=True)

scatter_ani.save('../Data/Evolution2D.gif', writer='pillow', fps=30)
plt.show()