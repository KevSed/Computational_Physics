import numpy as np
from matplotlib import pyplot as plt
from matplotlib import lines
linestyles = ['-', '--', '-.', ':']

plt.rc("text", usetex=True)

j,psi= np.loadtxt("Anfangszustand.txt", unpack=True)

fig = plt.figure(figsize=(5,5))
plt.plot(j, psi, linewidth=1)

#plt.ylim(-5,5)
#plt.xlim(-5, 5)
plt.ylabel(r"$\left|\Psi\left(\varepsilon_j,\tau=0\right)\right|^2$")
plt.xlabel(r"$j$")
plt.grid()

plt.savefig("Anfangszustand.pdf", bbox_inches="tight")
