from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import pandas as pd
from scipy.interpolate import griddata


fig = plt.figure()
ax = fig.gca(projection='3d')

# import data.
x, y, v = np.genfromtxt('Potential_V.txt', unpack='True')
xyv = {'x': x, 'y': y, 'v': v}
df = pd.DataFrame(xyv, index=range(len(xyv['x'])))
x1 = np.linspace(df['x'].min(), df['x'].max(), len(df['x'].unique()))
y1 = np.linspace(df['y'].min(), df['y'].max(), len(df['y'].unique()))
x2, y2 = np.meshgrid(x1, y1)
z2 = griddata((df['x'], df['y']), df['v'], (x2, y2), method='cubic')
# Plot the surface.
surf = ax.plot_surface(x2, y2, z2, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

# Customize the z axis.
#ax.set_zlim(-1.01, 1.01)
#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Potential V(x,y)')
plt.title('Three magnets system')

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.savefig('Potential.pdf')
plt.show()
