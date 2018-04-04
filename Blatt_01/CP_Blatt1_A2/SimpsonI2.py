import numpy as np
from matplotlib import pyplot as plt
from matplotlib import lines
linestyles = ['-', '--', '-.', ':']

plt.rc("text", usetex=True)

N,I= np.loadtxt("I2_Simpson.txt", unpack=True)

fig = plt.figure(figsize=(6,4))
plt.plot(N, I)

plt.ylim(0.28, 0.45)
plt.xlim(0, 2048)
plt.ylabel(r"$I_2$")
plt.xlabel(r"$N$")
plt.grid()




plt.savefig("I2Simpson.pdf", bbox_inches="tight")
