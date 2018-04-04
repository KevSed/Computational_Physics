import numpy as np
from matplotlib import pyplot as plt
from matplotlib import lines
linestyles = ['-', '--', '-.', ':']

plt.rc("text", usetex=True)

x, E= np.loadtxt("Energie.txt", unpack=True)

fig = plt.figure(figsize=(6,4))
plt.plot(x, E, label=r"$E(x)$")
plt.ylim(-3, -1)
plt.xlim(-3.1, 3.1)
plt.ylabel(r"$E\left[\frac{q_1 q_2}{4\pi\varepsilon_0}\right]$")
plt.xlabel(r"$x[a]$")
plt.grid()
plt.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2)

plt.savefig("Energie.pdf", bbox_inches="tight")
