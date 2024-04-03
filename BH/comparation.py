import matplotlib.pyplot as plt
import numpy as np

# Datos
N_bodies = [2, 10, 100, 1000, 10000]
time_barneshut = [0.042, 0.11, 2.8, 136.2, 4282.9]
time_Nbody = [0.06, 0.22, 7.54, 713.9, 70817.8]

# Complejidades esperadas

# Crear la gráfica
plt.figure(figsize=(10, 6))
plt.plot(N_bodies, time_barneshut, marker='.', label='Barnes-Hut', color='r', linewidth=1)
plt.plot(N_bodies, time_Nbody, marker='.', label='N-body', color='b', linewidth=1)

# Etiquetas y título
plt.title('Strong Scaling Comparison: Barnes-Hut vs Brute Force')
plt.xlabel('Number of Bodies [N]')
plt.ylabel('Time [s]')
plt.xscale('log')
plt.yscale('log')
plt.xticks(N_bodies, N_bodies)
plt.legend()

# Mostrar la gráfica
plt.tight_layout()
plt.show()
