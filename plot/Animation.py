from numpy import loadtxt
import matplotlib.pyplot as plt
plt.style.use('dark_background')
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
from sys import argv

N = int(argv[1])    # Number of particles
dt = float(argv[2])    # Size step
jump = int(argv[3])    # Size jump

x, y, z, m = loadtxt("../Data/Evolution.txt", unpack=1)    # Evolution
steps = x.size // N # Printed steps

# Extract the initial positions of the particles
initial_x = x[:N]
initial_y = y[:N]
initial_z = z[:N]

def update_scatter(num, scatter, message):
    message.set_text(f"Time = {num*jump*dt: .2f} years")
    scatter._offsets3d = (x[num*N:(num+1)*N], y[num*N:(num+1)*N], z[num*N:(num+1)*N])
    return scatter, message

# Attaching 3D axis to the figure
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

message = ax.text2D(0.00, 0.9, "", transform=ax.transAxes)

# Creating a scatter plot for N bodies
scatter = ax.scatter(initial_x, initial_y, initial_z, c='white', marker='.')

boundaries = int(argv[4])

# Setting the axes properties
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

scatter_ani.save('Evolution.gif', writer='pillow', fps=30)
plt.show()