from numpy import loadtxt
import matplotlib.pyplot as plt
plt.style.use('dark_background')
import mpl_toolkits.mplot3d.axes3d as p3
from sys import argv

N = int(argv[1])        # Number of particles

x, y, z, m = loadtxt("../Data/Evolution.txt", unpack=1, max_rows=N)
steps = x.size // N     # Printed steps

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Creating a scatter plot for N bodies
scatter = ax.scatter(x, y, z, c='white', marker='.')

boundaries = int(argv[2])

#Setting the axes properties
if boundaries:
    ax.set_xlim3d([-boundaries, boundaries])
    ax.set_ylim3d([-boundaries, boundaries])
    ax.set_zlim3d([-boundaries, boundaries])
    ax.set_xlabel('x [au]')
    ax.set_ylabel('y [au]')
    ax.set_zlabel('z [au]')

plt.grid(False)
plt.axis('off')
plt.show()
