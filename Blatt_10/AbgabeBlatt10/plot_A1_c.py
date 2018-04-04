import numpy as np
from matplotlib import pyplot as plt
from matplotlib import lines
linestyles = ['-', '--', '-.', ':']

plt.rc("text", usetex=True)

x,px= np.loadtxt("p.txt", unpack=True)
a= np.loadtxt("bi_ruck.txt", unpack=True)
y= np.loadtxt("bii_ruck.txt", unpack=True)
z= np.loadtxt("mt_ruck.txt", unpack=True)



fig = plt.figure(figsize=(7,7))

plt.plot(x,px) 
#plt.hist(a,bins=40)
plt.hist(a,bins=40, normed="true")
plt.xlabel("x")
plt.ylabel("p(x)")
plt.savefig("p.pdf")
plt.close();


plt.plot(x,px)
#plt.hist(y,bins=40) 
plt.hist(y,bins=40, normed="true")
plt.xlabel("x")
plt.ylabel("p(x)")
plt.savefig("p_ii.pdf")
plt.close();


plt.plot(x,px)
#plt.hist(z,bins=40) 
plt.hist(z,bins=40, normed="true")
plt.xlabel("x")
plt.ylabel("p(x)")
plt.savefig("p_mt.pdf")
plt.close();
