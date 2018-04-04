import numpy as np
from matplotlib import pyplot as plt
from matplotlib import lines
from mpl_toolkits import mplot3d
linestyles = ['-', '--', '-.', ':']

plt.rc("text", usetex=True)

x= np.loadtxt("bi.txt", unpack=True)
y= np.loadtxt("bii.txt", unpack=True)
z= np.loadtxt("mt.txt", unpack=True)

xt1,xt2= np.loadtxt("bi_test1.txt", unpack=True)
yt1,yt2= np.loadtxt("bii_test1.txt", unpack=True)
zt1,zt2= np.loadtxt("mt_test.txt", unpack=True)

at1,at2,at3= np.loadtxt("bi_test2.txt", unpack=True)
bt1,bt2,bt3= np.loadtxt("bii_test2.txt", unpack=True)
ct1,ct2,ct3= np.loadtxt("mt_test2.txt", unpack=True)

fig = plt.figure(figsize=(7,7))

plt.hist(x, bins=20, normed="true") 
plt.xlabel("Zufallszahl")
plt.ylabel("Haeufigkeit")
plt.savefig("bi_hist.pdf")
plt.close();

fig = plt.figure(figsize=(7,7))
plt.hist(y, bins=20, normed="true") 
plt.xlabel("Zufallszahl")
plt.ylabel("Haeufigkeit")
plt.savefig("bii_hist.pdf")
plt.close();

fig = plt.figure(figsize=(7,5))
plt.hist(z, bins=20, normed="true") 
plt.xlabel("Zufallszahl")
plt.ylabel("Haeufigkeit")
plt.savefig("mersenne_hist.pdf")
plt.close();

plt.scatter(xt1, xt2) 
plt.xlabel(r"$r_n$")
plt.ylabel(r"$r_{n+1}$")
plt.ylim(0,1)
plt.xlim(0,1)
plt.savefig("bi_test.pdf")
plt.close();

#3d
fig = plt.figure(figsize=(6,6))
ax = plt.axes(projection='3d')
ax.scatter3D(at1, at2, at3)

#plt.ylim(-3.1,3.1)
#plt.xlim(-3, 3)
ax.set_xlabel(r"$r_n$")
ax.set_ylabel(r"$r_{n+1}$")
ax.set_zlabel(r"$r_{n+2}$")
plt.savefig("bi_test2.pdf")
plt.grid()
for angle in range(0, 360):
    ax.view_init(30, angle)
    plt.draw()
    plt.pause(.001)
plt.close();

#3d
fig = plt.figure(figsize=(6,6))
ax = plt.axes(projection='3d')
ax.scatter3D(bt1, bt2, bt3)

#plt.ylim(-3.1,3.1)
#plt.xlim(-3, 3)
ax.set_xlabel(r"$r_n$")
ax.set_ylabel(r"$r_{n+1}$")
ax.set_zlabel(r"$r_{n+2}$")
plt.savefig("bii_test2.pdf")
plt.grid()
for angle in range(0, 360):
    ax.view_init(30, angle)
    plt.draw()
    plt.pause(.001)
plt.close();

#3d
fig = plt.figure(figsize=(6,6))
ax = plt.axes(projection='3d')
ax.scatter3D(ct1, ct2, ct3)

#plt.ylim(-3.1,3.1)
#plt.xlim(-3, 3)
ax.set_xlabel(r"$r_n$")
ax.set_ylabel(r"$r_{n+1}$")
ax.set_zlabel(r"$r_{n+2}$")
plt.savefig("mt_test2.pdf")
plt.grid()
for angle in range(0, 360):
    ax.view_init(30, angle)
    plt.draw()
    plt.pause(.001)
plt.close();

plt.scatter(yt1, yt2, marker='x', s=1) 
plt.xlabel(r"$r_n$")
plt.ylabel(r"$r_{n+1}$")
plt.ylim(0,1)
plt.xlim(0,1)
plt.savefig("bii_test.pdf")
plt.close();

plt.scatter(zt1, zt2, marker='x', s=1) 
plt.xlabel(r"$r_n$")
plt.ylabel(r"$r_{n+1}$")
plt.ylim(0,1)
plt.xlim(0,1)
plt.savefig("marsenne_test.pdf")
plt.close();
