import numpy as np
from matplotlib import pyplot as plt
from matplotlib import lines
linestyles = ['-', '--', '-.', ':']

plt.rc("text", usetex=True)

x, Fnum03, Fan, _= np.loadtxt("Kraft_03a.txt", unpack=True)
_, FnumE4, _, _= np.loadtxt("Kraft_E4a.txt", unpack=True)
_, FnumE15, _, _= np.loadtxt("Kraft_E15a.txt", unpack=True)

fig = plt.figure(figsize=(6,4))
plt.plot(x, FnumE15, 'c', linestyle='--' , dashes=(5,2), label=r"$F_{num}(x)$ mit $h=10^{-15}a$")
plt.plot(x, Fan, 'b', label=r"$F_{an}(x)$")
plt.plot(x, Fnum03, 'r', label=r"$F_{num}(x)$ mit $h=0.3a$")
plt.plot(x, FnumE4, 'g', linestyle='--' , dashes=(6,3), label=r"$F_{num}(x)$ mit $h=10^{-4}a$")

plt.ylim(-1.3, 1.3)
plt.xlim(-3.1, 3.1)
plt.ylabel(r"$F\left[\frac{q_1 q_2}{4\pi\varepsilon_0}\right]$")
plt.xlabel(r"$x[a]$")
plt.grid()
plt.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2)




plt.savefig("Kraft.pdf", bbox_inches="tight")
