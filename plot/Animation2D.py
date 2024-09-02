from numpy import loadtxt, array, fromfile, float64
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
    data = loadtxt("../input", max_rows=6)
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

# Data for xy
xi1 = 0
yi1 = 1
initial_x1, initial_y1 = read_data("../Data/Ev_0", N, [xi1, yi1])

# Data for xz
xi2 = 0
yi2 = 2
initial_x2, initial_y2 = read_data("../Data/Ev_0", N, [xi2, yi2])

def update_scatter(num, scatter1, message1, scatter2, message2):
    message1.set_text(f" t = {(num+1)*jump*dt: .2f} yr")
    x1, y1 = read_data("../Data/Ev_"+str(num), N, [xi1, yi1])
    scatter1.set_offsets(array([x1, y1]).transpose())
    
    x2, y2 = read_data("../Data/Ev_"+str(num), N, [xi2, yi2])
    scatter2.set_offsets(array([x2, y2]).transpose())
    
    return scatter1, message1, scatter2, message2

# Creating the figure with two subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

# Plot 1
message1 = ax1.text(0.0, 0.9, "", transform=ax1.transAxes)
scatter1 = ax1.scatter(initial_x1, initial_y1, c='white', marker='.', s=1)
boundaries = 7.0 * R

if boundaries:
    ax1.set_xlim([-boundaries, boundaries])
    ax1.set_ylim([-boundaries, boundaries])
    ax1.set_xlabel('x [au]')
    ax1.set_ylabel('y [au]')
    ax1.set_title('X-Y Plane')

# Plot 2
message2 = ax2.text(0.0, 0.9, "", transform=ax2.transAxes)
scatter2 = ax2.scatter(initial_x2, initial_y2, c='white', marker='.', s=1)

if boundaries:
    ax2.set_xlim([-boundaries, boundaries])
    ax2.set_ylim([-boundaries, boundaries])
    ax2.set_xlabel('x [au]')
    ax2.set_ylabel('y [au]')
    ax2.set_title('X-Z Plane')

# Creating the Animation object
scatter_ani = animation.FuncAnimation(fig, update_scatter, frames=steps, fargs=(scatter1, message1, scatter2, message2), interval=1, blit=True, repeat=False)

scatter_ani.save('../Data/Evolution2D.gif', writer='pillow', fps=10)
plt.tight_layout()
plt.show()
