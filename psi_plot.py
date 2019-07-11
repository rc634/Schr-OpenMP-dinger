import matplotlib.pyplot as plt
import numpy as np

max_plot_value = 1.0


NX1 = int(np.loadtxt("param.par")[0] + 2*np.loadtxt("param.par")[2])
NX2 = int(np.loadtxt("param.par")[1] + 2*np.loadtxt("param.par")[2])

init_re = np.fromfile('init_re.dat', dtype=np.float64).reshape((NX1,NX2))
init_im = np.fromfile('init_im.dat', dtype=np.float64).reshape((NX1,NX2))

fig, axs = plt.subplots(2,1,figsize=(8,10))
ax = axs[0]
re = ax.pcolormesh(init_re, cmap='bwr', vmin=-max_plot_value, vmax=max_plot_value)
fig.colorbar(re)
ax.set_aspect('equal')

ax = axs[1]
im = ax.pcolormesh(init_im, cmap='bwr', vmin=-max_plot_value, vmax=max_plot_value)
# fig.colorbar(im)
ax.set_aspect('equal')
plt.title("Real Value Top, Im Value Bottom")

#plt.show()

final_re = np.fromfile('final_re.dat', dtype=np.float64).reshape((NX1,NX2))
final_im = np.fromfile('final_im.dat', dtype=np.float64).reshape((NX1,NX2))

fig, axs = plt.subplots(2,1,figsize=(8,10))
ax = axs[0]
re = ax.pcolormesh(final_re, cmap='bwr', vmin=-max_plot_value, vmax=max_plot_value)
fig.colorbar(re)
ax.set_aspect('equal')

ax = axs[1]
im = ax.pcolormesh(final_im, cmap='bwr', vmin=-max_plot_value, vmax=max_plot_value)
# fig.colorbar(im)
ax.set_aspect('equal')
plt.title("Real Value Top, Im Value Bottom")

plt.show()
