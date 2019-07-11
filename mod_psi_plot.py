import matplotlib.pyplot as plt
import numpy as np

max_plot_value = 1.0


NX1 = int(np.loadtxt("param.par")[0] + 2*np.loadtxt("param.par")[2])
NX2 = int(np.loadtxt("param.par")[1] + 2*np.loadtxt("param.par")[2])

init = np.fromfile('init_modsqr.dat', dtype=np.float64).reshape((NX1,NX2))

fig, axs = plt.subplots(1,1,figsize=(8,10))
ax = axs
re = ax.pcolormesh(init, cmap='rainbow', vmin=0, vmax=max_plot_value)
fig.colorbar(re)
ax.set_aspect('equal')

# fig.colorbar(im)
plt.title(r"Initial probablility density $|\psi(0)|^2$")

#plt.show()

final = np.fromfile('final_modsqr.dat', dtype=np.float64).reshape((NX1,NX2))

fig, axs = plt.subplots(1,1,figsize=(8,10))
ax = axs
re = ax.pcolormesh(final, cmap='rainbow', vmin=0, vmax=max_plot_value)
fig.colorbar(re)
ax.set_aspect('equal')

# fig.colorbar(im)
plt.title(r"Final $|\psi(t_0)|^2$")

plt.show()
