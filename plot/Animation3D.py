from numpy import loadtxt
import matplotlib.pyplot as plt
plt.style.use('dark_background')
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation

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


N, dt, jump= read_parameters([0, 3, 6])
N = int(N)
jump = int(jump)

x, y, z, m = loadtxt("../Data/Evolution.txt", unpack=1)    # Evolution
steps = x.size // N # Printed steps

# Extract the initial positions of the particles
initial_x = x[:N]
initial_y = y[:N]
initial_z = z[:N]

def update_scatter(num, scatter, message):
    message.set_text(f"t = {num*jump*dt: .2f} yr")
    scatter._offsets3d = (x[num*N:(num+1)*N], y[num*N:(num+1)*N], z[num*N:(num+1)*N])
    return scatter, message


# Attaching 3D axis to the figure
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

message = ax.text2D(0.00, 0.9, "", transform=ax.transAxes)

# Creating a scatter plot for N bodies
scatter = ax.scatter(initial_x, initial_y, initial_z, c='white', marker='.', s=1)

# Setting the axes properties
boundaries = 1.0
if boundaries:
    ax.set_xlim3d([-boundaries, boundaries])
    ax.set_ylim3d([-boundaries, boundaries])
    ax.set_zlim3d([-boundaries, boundaries])
    ax.set_xlabel('x [au]')
    ax.set_ylabel('y [au]')
    ax.set_zlabel('z [au]')

# Creating the Animation object
scatter_ani = animation.FuncAnimation(fig, update_scatter, frames=steps, fargs=(scatter, message),
                                      interval=1, blit=True)

scatter_ani.save('../Data/Evolution.gif', writer='pillow', fps=10)
plt.show()