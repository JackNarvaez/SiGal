from numpy import loadtxt, array, fromfile, float64, int32
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
    data = fromfile(file_path, dtype=float64, count=3 * N)
    ids = fromfile(file_path, dtype=int32, offset=8 * 7 * N)
    x = data[coords[0]:3*N+coords[0]:3]
    y = data[coords[1]:3*N+coords[1]:3]
    
    return x, y, ids

warm_colors = ['#FFFFFF', '#FFA500', '#FF4500', '#FFD700', '#FF1493']
cool_colors = ['#FFFFFF', '#00FFFF', '#4169E1', '#800080', '#00FF00']
color_map   = warm_colors

# --------------------- Parameters --------------------- #

N, R, dt, steps, jump = read_parameters([0, 2, 3, 4, 5])
N = int(N)
jump = int(jump)
steps = int(steps)//jump
Dark_Matter = 2   # <--- MODIFY THIS ACCORDING TO THE SYSTEM!!!
# ------------------------------------------------------ #

# Data for xy
xi1 = 0
yi1 = 1
initial_x1, initial_y1, id_1= read_data("../Data/Ev_0", N, [xi1, yi1])

# Data for xz
xi2 = 0
yi2 = 2
initial_x2, initial_y2, id_2 = read_data("../Data/Ev_0", N, [xi2, yi2])

def update_scatter(num, scatter1, message1, scatter2, message2, exclude_id):
    message1.set_text(f" t = {(num+1)*jump*dt: .2f} yr")
    x1, y1, id1 = read_data("../Data/Ev_" + str(num), N, [xi1, yi1])
    mask1 = id1 != exclude_id
    x1_filtered = x1[mask1]
    y1_filtered = y1[mask1]
    id1_filtered = id1[mask1]
    colors1 = [color_map[i % len(color_map)] for i in id1_filtered]
    scatter1.set_offsets(array([x1_filtered, y1_filtered]).transpose())
    scatter1.set_color(colors1)

    x2, y2, id2 = read_data("../Data/Ev_" + str(num), N, [xi2, yi2])
    mask2 = id2 != exclude_id
    x2_filtered = x2[mask2]
    y2_filtered = y2[mask2]
    id2_filtered = id2[mask2]
    colors2 = [color_map[i % len(color_map)] for i in id2_filtered]
    scatter2.set_offsets(array([x2_filtered, y2_filtered]).transpose())
    scatter2.set_color(colors2)

    return scatter1, message1, scatter2, message2


# Creating the figure with two subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

# Plot 1
mask1 = id_1 != Dark_Matter
initial_x1_filtered = initial_x1[mask1]
initial_y1_filtered = initial_y1[mask1]
id_1_filtered = id_1[mask1]
colors1 = [color_map[i % len(color_map)] for i in id_1_filtered]

message1 = ax1.text(0.0, 0.9, "", transform=ax1.transAxes)
scatter1 = ax1.scatter(initial_x1_filtered, initial_y1_filtered, c=colors1, marker='.', s=1)

boundaries = 2.0 * R

if boundaries:
    ax1.set_xlim([-boundaries, boundaries])
    ax1.set_ylim([-boundaries, boundaries])
    ax1.set_xlabel('x [au]')
    ax1.set_ylabel('y [au]')
    ax1.set_title('X-Y Plane')

# Plot 2

mask2 = id_2 != Dark_Matter
initial_x2_filtered = initial_x2[mask2]
initial_y2_filtered = initial_y2[mask2]
id_2_filtered = id_2[mask2]
colors2 = [color_map[i % len(color_map)] for i in id_2_filtered]

message2 = ax2.text(0.0, 0.9, "", transform=ax2.transAxes)
scatter2 = ax2.scatter(initial_x2_filtered, initial_y2_filtered, c=colors2, marker='.', s=1)

if boundaries:
    ax2.set_xlim([-boundaries, boundaries])
    ax2.set_ylim([-boundaries, boundaries])
    ax2.set_xlabel('x [au]')
    ax2.set_ylabel('y [au]')
    ax2.set_title('X-Z Plane')

# Creating the Animation object
scatter_ani = animation.FuncAnimation(fig, update_scatter, frames=steps, fargs=(scatter1, message1, scatter2, message2, Dark_Matter), interval=1, blit=True, repeat=False)

scatter_ani.save('../Data/Evolution2D.gif', writer='pillow', fps=10)
plt.tight_layout()
plt.show()