import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D  # Importa Axes3D para gráficos 3D
import pandas as pd

num_steps = 1200  # Ajusta esto al número de pasos de tu simulación

# Función para leer las posiciones de las partículas de un archivo CSV
def load_particle_positions(step):
    filename = f'particulas_posiciones_{step}.csv'
    data = pd.read_csv(filename)
    x_positions = data['X'].values
    y_positions = data['Y'].values
    z_positions = data['Z'].values
    return x_positions, y_positions, z_positions

# Configuración inicial de la figura y los ejes de matplotlib para la animación en 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')  # Configura el subplot para 3D

# Configura los límites de los ejes basado en tu simulación
ax.set_xlim([-100, 100])
ax.set_ylim([-100, 100])
ax.set_zlim([-100, 100])

# Inicializa un gráfico de puntos para las partículas en 3D
particles, = ax.plot([], [], [], 'bo', ms=1)  # 'bo' crea puntos azules, ms es el tamaño del punto

# Función de inicialización para FuncAnimation
def init():
    particles.set_data([], [])
    particles.set_3d_properties([])  # Inicializa las propiedades 3D
    return particles,

# Función de animación que se llama secuencialmente
def animate(step):
    x, y, z = load_particle_positions(step)
    particles.set_data(x, y)
    particles.set_3d_properties(z)  # Actualiza las propiedades 3D
    return particles,

# Crea la animación usando FuncAnimation
ani = animation.FuncAnimation(fig, animate, init_func=init, frames=num_steps, interval=20, blit=False)

plt.show()

# Para guardar la animación, puedes descomentar la siguiente línea:
ani.save('simul_3D.mp4', writer='ffmpeg', fps=100)

