import numpy as np
from matplotlib import pyplot as plt

plt.rc("text", usetex=True)

x, _, _, Frel03= np.loadtxt("Kraft_03a.txt", unpack=True)
_, _, _, Frel4= np.loadtxt("Kraft_E4a.txt", unpack=True)
_, _, _, Frel15= np.loadtxt("Kraft_E15a.txt", unpack=True)

fig = plt.figure(figsize=(5,10))

fig.add_subplot(3, 1, 1)

plt.plot(x, Frel03)
plt.xlim(-1, 1)
plt.ylim(-10, 10)
plt.ylabel(r"$\frac{F_{num}(x)-F_{an}(x)}{F_an(x)}$ mit $h=0.3 a$")
plt.grid()

fig.add_subplot(3, 1, 2)

plt.plot(x, Frel4)
plt.xlim(-1, 1)
plt.ylim(-10e-6, 10e-6)
plt.ylabel(r"$\frac{F_{num}(x)-F_{an}(x)}{F_an(x)}$ mit $h=10^{-4} a$")
plt.grid()

fig.add_subplot(3, 1, 3)

plt.plot(x, Frel15)
plt.xlim(-1, 1)
plt.ylim(-10, 10)
plt.ylabel(r"$\frac{F_{num}(x)-F_{an}(x)}{F_an(x)}$ mit $h=10^{-15} a$")
plt.grid()


plt.xlabel(r"$x[a]$")


plt.savefig("relFehler_zoom.pdf", bbox_inches="tight")
