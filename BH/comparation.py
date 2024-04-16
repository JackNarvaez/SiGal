import matplotlib.pyplot as plt
import numpy as np

# Datos existentes
N_bodies = np.array([10, 100, 1000, 10000])
time_barneshut = np.array([0.11, 2.8, 136.2, 4282.9])
time_Nbody = np.array([0.22, 7.54, 713.9, 70817.8])

# Cálculos para ajustar las curvas teóricas
k_Nbody = np.mean(time_Nbody / (N_bodies**2))
k_BarnesHut = np.mean(time_barneshut / (N_bodies * np.log(N_bodies)))

# Generar datos teóricos para las curvas
N_teoretical = np.linspace(min(N_bodies), max(N_bodies), 100)
time_teoretical_N2 = k_Nbody * N_teoretical**2
time_teoretical_NlogN = k_BarnesHut * N_teoretical * np.log(N_teoretical)

# Crear la gráfica con las curvas teóricas añadidas
plt.figure(figsize=(8, 6))
plt.plot(N_bodies, time_barneshut, marker='.', label='Barnes-Hut', color='r', linewidth=1)
plt.plot(N_bodies, time_Nbody, marker='.', label='N-body (Brute Force)', color='b', linewidth=1)

# Etiquetas y título
plt.title('Strong Scaling Comparison', fontsize=20)
plt.xlabel('N', fontsize=20)
plt.ylabel('T [s]', fontsize=20)
plt.xscale('log')
plt.yscale('log')
plt.xticks(N_bodies, labels=N_bodies, fontsize=12)
plt.yticks(fontsize=20)
plt.legend(fontsize=20)

# Ajustes adicionales para mejorar la visualización
plt.tight_layout()
# Mostrar la gráfica
plt.show()
