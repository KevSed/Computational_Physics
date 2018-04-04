import numpy as np
from matplotlib import pyplot as plt
from matplotlib import lines
linestyles = ['-', '--', '-.', ':']

plt.rc("text", usetex=True)

j,psi= np.loadtxt("V=0", unpack=True)
j10,psi10= np.loadtxt("V=10", unpack=True)
j30,psi30= np.loadtxt("V=30", unpack=True)
j50,psi50= np.loadtxt("V=50", unpack=True)

fig = plt.figure(figsize=(10,5))
plt.plot(j, psi, linewidth=1, label="V=0")
plt.plot(j10, psi10, linewidth=1, label="V=10")
plt.plot(j30, psi30, linewidth=1, label="V=30")
plt.plot(j50, psi50, linewidth=1, label="V=50")

#plt.ylim(-5,5)
#plt.xlim(-5, 5)
plt.ylabel(r"$\left|\Psi\left(\varepsilon_j,\tau=1\right)\right|^2$")
plt.xlabel(r"$\varepsilon_j$")
plt.grid()
plt.legend(loc="best")


plt.savefig("Zeitentwicklung.pdf", bbox_inches="tight")
