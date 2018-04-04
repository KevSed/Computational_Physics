import numpy as np
from matplotlib import pyplot as plt
from matplotlib import lines
linestyles = ['-', '--', '-.', ':']

plt.rc("text", usetex=True)

N,I= np.loadtxt("I1_Mittelpunkt.txt", unpack=True)

fig = plt.figure(figsize=(6,4))
plt.plot(N, I)

plt.ylim(0, 0.23)
plt.xlim(0, 4096)
plt.ylabel(r"$I_1$")
plt.xlabel(r"$N$")
plt.grid()




plt.savefig("I1Mittel.pdf", bbox_inches="tight")
