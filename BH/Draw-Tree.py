import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def draw_cube(ax, corners, line_width=0.5):
    min_x, min_y, min_z, max_x, max_y, max_z = corners
    for s, e in [((min_x, min_y, min_z), (max_x, min_y, min_z)),
                 ((min_x, min_y, min_z), (min_x, max_y, min_z)),
                 ((min_x, min_y, min_z), (min_x, min_y, max_z)),
                 ((max_x, max_y, max_z), (min_x, max_y, max_z)),
                 ((max_x, max_y, max_z), (max_x, min_y, max_z)),
                 ((max_x, max_y, max_z), (max_x, max_y, min_z)),
                 ((min_x, max_y, max_z), (min_x, max_y, min_z)),
                 ((min_x, max_y, max_z), (max_x, max_y, max_z)),
                 ((max_x, min_y, min_z), (max_x, max_y, min_z)),
                 ((max_x, min_y, min_z), (max_x, min_y, max_z)),
                 ((min_x, max_y, min_z), (min_x, min_y, min_z)),
                 ((max_x, max_y, min_z), (max_x, max_y, max_z))]:
        ax.plot3D(*zip(s, e), color="k", lw=line_width)

def draw_rectangle(ax, corners, line_width=0.5):
    min_x, min_y, min_z, max_x, max_y, max_z = corners
    xy = [(min_x, min_y), (max_x, min_y), (max_x, max_y), (min_x, max_y), (min_x, min_y)]
    ax.plot(*zip(*xy), color="k", lw=line_width)


bodies_df = pd.read_csv("positions_0.csv")

fig_3d = plt.figure(figsize=(10, 8))
ax3d = fig_3d.add_subplot(111, projection='3d')

with open("octants.txt", "r") as file:
    for line in file:
        corners = tuple(map(float, line.strip().split(',')))
        draw_cube(ax3d, corners, line_width=0.9)

ax3d.scatter(bodies_df['X'], bodies_df['Y'], bodies_df['Z'], color='r', s=5)

ax3d.set_title('3D View')
plt.savefig('3D-Tree.png',dpi=600)
plt.show()


#_________________________________________________________________#
#____________________________2D-PLOT___________________________________#
#__________________________________________________________________


fig_2d, axs_2d = plt.subplots(1, 3, figsize=(18, 6))
ax_xy, ax_xz, ax_yz = axs_2d

with open("octants.txt", "r") as file:
    for line in file:
        corners = tuple(map(float, line.strip().split(',')))
        draw_rectangle(ax_xy, corners, line_width=1)  
        draw_rectangle(ax_xz, (corners[0], corners[2], 0, corners[3], corners[5], 0), line_width=1)  
        draw_rectangle(ax_yz, (corners[1], corners[2], 0, corners[4], corners[5], 0), line_width=1)  

ax_xy.scatter(bodies_df['X'], bodies_df['Y'], color='r', s=5)
ax_xz.scatter(bodies_df['X'], bodies_df['Z'], color='r', s=5)
ax_yz.scatter(bodies_df['Y'], bodies_df['Z'], color='r', s=5)

ax_xy.set_title('XY Plane')
ax_xz.set_title('XZ Plane')
ax_yz.set_title('YZ Plane')

plt.tight_layout()
plt.savefig('2D-Tree.png',dpi=600)
plt.show()
