import numpy as np
from matplotlib import pyplot as plt
from matplotlib import lines
linestyles = ['-', '--', '-.', ':']

plt.rc("text", usetex=True)

t0,m0= np.loadtxt("Data_h0", unpack=True)
t1,m1= np.loadtxt("Data_h01", unpack=True)
t5,m5= np.loadtxt("Data_h05", unpack=True)

h05,mh05= np.loadtxt("Data_t05", unpack=True)
h1,mh1= np.loadtxt("Data_t10", unpack=True)
h15,mh5= np.loadtxt("Data_t15", unpack=True)

fig = plt.figure(figsize=(5,5))
plt.plot(t0, m0, "b.", markersize=0.5,label="h=0")
plt.plot(t1, m1, "r.", markersize=0.5,label="h=0.1")
plt.plot(t5, m5, "g.", markersize=0.5,label="h=0.5")

#plt.ylim(-5,5)
#plt.xlim(-5, 5)
plt.ylabel(r"$m$")
plt.xlabel(r"$t$")
plt.grid()
plt.legend(loc="best")
plt.savefig("mag_tvari.pdf", bbox_inches="tight")


fig = plt.figure(figsize=(5,5))
plt.plot(h05, mh05, "b.", markersize=0.5, label="t=0.5")
plt.plot(h1, mh1, "k.", markersize=0.5, label="t=1.0")
plt.plot(h15, mh5, "r.", markersize=0.5, label="t=1.5")

#plt.ylim(-5,5)
#plt.xlim(-5, 5)
plt.ylabel(r"$m$")
plt.xlabel(r"$h$")
plt.grid()
plt.legend(loc="best")
plt.savefig("mag_hvari.pdf", bbox_inches="tight")
