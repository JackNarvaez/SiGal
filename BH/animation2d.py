import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pandas as pd

num_steps = 1000

def load_particle_positions(step):
    filename = f'/home/yo/Documents/NBodySimulations/BH/files/positions_{step}.csv'
    data = pd.read_csv(filename)
    x_positions = data['X'].values
    y_positions = data['Y   '].values
    # No necesitamos z_positions para la animación 2D
    return x_positions, y_positions

fig = plt.figure()
ax = fig.add_subplot(111)  
ax.set_xlim([-20, 60])
ax.set_ylim([-60, 20])
ax.set_facecolor('black')  
particles, = ax.plot([], [], 'wo', ms=0.5)  
def init():
    particles.set_data([], [])
    return particles,

def animate(step):
    x, y = load_particle_positions(step)
    particles.set_data(x, y)  
    return particles,

ani = animation.FuncAnimation(fig, animate, init_func=init, frames=num_steps, interval=10, blit=False)

plt.show()

# Guardar como video 2D
ani.save('/home/yo/Documents/NBodySimulations/BH/files/simul_2D.gif', writer='pillow', fps=300) 
