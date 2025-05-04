import matplotlib.pyplot as plt
from Data   import read_parameters,read_full_data
from numpy  import sum,sqrt,zeros,linalg,newaxis,cross,arange

def kin_energ(m, vel):
    v2 = sum(vel**2, axis=1)
    return 0.5*sum(m * v2)

def pot_energ(m, pos):
    N = len(m)
    U = 0.0
    for i in range(N):
        r_ij = pos[i] - pos[i+1:]
        dist = linalg.norm(r_ij, axis=1)
        U -= sum(m[i] * m[i+1:] / dist)
    return U

def momentum(m, pos, vel):
    m = m[:, newaxis]
    P = sum(m*vel,axis=0)
    L = sum(m*cross(pos, vel), axis=0)
    return P, L

# --------------------- Parameters --------------------- #

N, R, dt, steps, jump = read_parameters([0, 2, 3, 4, 5])
N = int(N)
jump = int(jump)
steps = int(steps)//jump
# ------------------------------------------------------ #

time = dt*jump*arange(steps)
Ekin = zeros(steps)
Egrv = zeros(steps)
pmom = zeros((steps, 3))
lmom = zeros((steps, 3))

for ii in range(steps):
    pos, vel, mass, _ = read_full_data("../Data/Ev_"+str(ii), N)
    Ekin[ii] = kin_energ(mass, vel)
    Egrv[ii] = pot_energ(mass, pos)
    pmom[ii], lmom[ii] = momentum(mass, pos, vel)
    print(f'step:\t {ii}')

P = sqrt(sum((pmom-pmom[0])*(pmom-pmom[0]), axis=1))
L = sqrt(sum((lmom-lmom[0])*(lmom-lmom[0]), axis=1))

fig, ax = plt.subplots(figsize=(8,6))
ax.plot(time, Ekin,     ls="-", lw=3, c="forestgreen", label='Kinetic energy')
ax.plot(time, Egrv,     ls="-", lw=3, c="skyblue"    , label='Gravit. energy')
ax.plot(time, Ekin+Egrv,ls="-", lw=3, c="crimson"    , label='Total energy')
ax.set_xlim(time[0], time[-1])
ax.set_xlabel("time")
ax.legend()
plt.savefig("Energy.png")
plt.show()

fig, ax = plt.subplots(figsize=(8,6))
ax.plot(time, P,     ls="-", lw=3, c="forestgreen", label='Linear Mom.')
ax.plot(time, L,     ls="-", lw=3, c="crimson"    , label='Angular Mom.')
ax.set_xlim(time[0], time[-1])
ax.set_xlabel("time")
ax.legend()
plt.savefig("Momenta.png")
plt.show()