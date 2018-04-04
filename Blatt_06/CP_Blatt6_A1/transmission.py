import numpy as np
from matplotlib import pyplot as plt
from matplotlib import lines
linestyles = ['-', '--', '-.', ':']

plt.rc("text", usetex=True)

t,tr= np.loadtxt("V=0_trans", unpack=True)
t10,tr10= np.loadtxt("V=10_trans", unpack=True)
t30,tr30= np.loadtxt("V=30_trans", unpack=True)
t50,tr50= np.loadtxt("V=50_trans", unpack=True)

fig = plt.figure(figsize=(10,5))
plt.plot(t, tr, linewidth=1, label="V=0")
plt.plot(t10, tr10, linewidth=1, label="V=10")
plt.plot(t30, tr30, linewidth=1, label="V=30")
plt.plot(t50, tr50, linewidth=2, label="V=50")

plt.ylim(0,1.1)
plt.xlim(0, 1)
plt.ylabel(r"$Transmissionskoeffizient$")
plt.xlabel(r"$\tau$")
plt.grid()
plt.legend(loc="best")


plt.savefig("Transmission.pdf", bbox_inches="tight")
