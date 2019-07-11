import matplotlib.pyplot as plt
import numpy as np
import time

max_plot_value = 1.0


NX1 = int(np.loadtxt("param.par")[0] + 2*np.loadtxt("param.par")[2])
NX2 = int(np.loadtxt("param.par")[1] + 2*np.loadtxt("param.par")[2])

psi = []
FPS = 10.

i = 0
while True:
	try:
		if i<10:
			number = '0' + str(i)
		else:
			number = str(i)
		psi.append(np.fromfile('modsqr_'+number+'.dat', dtype=np.float64).reshape((NX1,NX2)))
		i += 1
	except:
		break

fig, ax = plt.subplots(1,1,figsize=(8,10))
ax.set_aspect('equal')

i = 0
while True:
	try:
		if i<10:
			number = '0' + str(i)
		else:
			number = str(i)
		i = int(i)
		re = ax.pcolormesh(psi[i], cmap='rainbow', vmin=0, vmax=max_plot_value)
		if i==0:
		    fig.colorbar(re)
		#plt.pcolormesh(psi[i], cmap='rainbow', vmin=0, vmax=max_plot_value)
		plt.title(r"Final $|\psi(t_{"+number+r"})|^2$")
		plt.savefig("movie/gaussian_" + number + ".png")
		plt.pause(1./FPS)
		i += 1
	except:
		break





