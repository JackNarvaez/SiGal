from numpy import loadtxt, array, fromfile, float64, array
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
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
    z = data[coords[2]:3*N+coords[2]:3]
    return x, y, z

N, R, dt, steps, jump = read_parameters([0, 2, 3, 4, 5])
N = int(N)
jump = int(jump)
steps = int(steps)//jump

#plot
xi = 0
yi = 1
zi = 2

# Extract the initial positions of the particles
initial_x, initial_y, initial_z = read_data("../Data/Ev_0", N, [xi,yi,zi])

def update_scatter(num, scatter, message):
    message.set_text(f"t = {(num+1)*jump*dt: .2f} yr")
    x, y, z = read_data("../Data/Ev_"+str(num), N, [xi,yi,zi])
    scatter._offsets3d = array([x, y, z])
    return scatter, message

# Attaching 3D axis to the figure
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

message = ax.text2D(0.00, 0.9, "", transform=ax.transAxes)

# Creating a scatter plot for N bodies
scatter = ax.scatter(initial_x, initial_y, initial_z, c='white', marker='.', s=1)

# Setting the axes properties
boundaries = 5.0*R
if boundaries:
    ax.set_xlim3d([-boundaries, boundaries])
    ax.set_ylim3d([-boundaries, boundaries])
    ax.set_zlim3d([-boundaries, boundaries])
    ax.set_xlabel('x [au]')
    ax.set_ylabel('y [au]')
    ax.set_zlabel('z [au]')

# Creating the Animation object
scatter_ani = animation.FuncAnimation(fig, update_scatter, frames=steps, fargs=(scatter, message), interval=1, blit=True, repeat=False)

scatter_ani.save('../Data/Evolution3D.gif', writer='pillow', fps=10)
plt.show()