import numpy as np
from matplotlib import pyplot as plt
from matplotlib import lines
linestyles = ['-', '--', '-.', ':']

plt.rc("text", usetex=True)

r,x= np.loadtxt("logistische_Abb", unpack=True)
rk,xk= np.loadtxt("kubische_Abb", unpack=True)

fig = plt.figure(figsize=(5,5))
plt.plot(r, x, "b.", markersize=0.1)

#plt.ylim(-5,5)
#plt.xlim(-5, 5)
plt.ylabel(r"$\left|\Psi\left(\varepsilon_j,\tau=0\right)\right|^2$")
plt.xlabel(r"$j$")
plt.grid()

plt.savefig("log_Abb.pdf", bbox_inches="tight")


fig = plt.figure(figsize=(5,5))
plt.plot(rk, xk, "b.", markersize=0.1)

#plt.ylim(-5,5)
#plt.xlim(-5, 5)
plt.ylabel(r"$\left|\Psi\left(\varepsilon_j,\tau=0\right)\right|^2$")
plt.xlabel(r"$j$")
plt.grid()

plt.savefig("kub_Abb.pdf", bbox_inches="tight")