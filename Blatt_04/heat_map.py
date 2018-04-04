import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 14})
mpl.rc('text', usetex=True)

#V1 = np.genfromtxt('Heat.txt')[::-1]
V2 = np.genfromtxt('Heat_b25.txt')[::-1]
'''
N = len(V1)

plt.figure(figsize=(7,7))
plt.title(r'$\beta=0.1$')
#plt.title('Potential $V(x,y)$')
plt.imshow(V1, interpolation='kaiser') #, cmap=mpl.cm.get_cmap('jet')
#plt.plot(0, 1, 'bo')
#plt.plot(-np.sqrt(3)/2,-0.5,'bo', label='Magnet')
#plt.plot( np.sqrt(3)/2,-0.5,'bo')
plt.xlabel(r'$y$')
plt.ylabel(r'$x$')
plt.tight_layout()
plt.gca().set_xticks(np.arange(-1 , len(V1), len(V1)/4));
plt.gca().set_yticks(np.arange(len(V1)-1, -1, -len(V1)/4));
plt.gca().set_xticklabels(np.arange(-2, 3, 1));
plt.gca().set_yticklabels(np.arange(-2, 3, 1));
#plt.show()
plt.savefig('heat_map_b1.pdf')

plt.cla()
plt.clf()
'''

N = len(V2)

plt.figure(figsize=(7,7))
plt.title(r'$\beta=0.25$')
#plt.title('Potential $V(x,y)$')
plt.imshow(V2, interpolation='kaiser') #, cmap=mpl.cm.get_cmap('jet')
#plt.plot(0, 1, 'bo')
#plt.plot(-np.sqrt(3)/2,-0.5,'bo', label='Magnet')
#plt.plot( np.sqrt(3)/2,-0.5,'bo')
plt.xlabel(r'$y$')
plt.ylabel(r'$x$')
plt.tight_layout()
plt.gca().set_xticks(np.arange(-1 , len(V2), len(V2)/4));
plt.gca().set_yticks(np.arange(len(V2)-1, -1, -len(V2)/4));
plt.gca().set_xticklabels(np.arange(-2, 3, 1));
plt.gca().set_yticklabels(np.arange(-2, 3, 1));
plt.colorbar()
#plt.show()
plt.savefig('heat_map_b25.pdf')
