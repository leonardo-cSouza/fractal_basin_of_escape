import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.gridspec as gridspec

plt.rc('text', usetex = True)
plt.rc('font', **{'family' : "sans-serif"})
params = {'text.latex.preamble' : [r'\usepackage{siunitx}', r'\usepackage{amsmath}']}
plt.rcParams.update(params)
plt.rcParams['axes.formatter.use_locale'] = True

# import the data
DATA = np.loadtxt("exits.dat") 

x = DATA[:,0]; y = DATA[:,1]; z = DATA[:,2]; t = DATA[:,3]; 

S = 1.0

fig = plt.figure(figsize=(8, 7))

cmap = colors.ListedColormap(['lime', 'red'])
#bounds take the third column of the file, here is z
bounds=[0, np.pi,2.0*np.pi]
norm = colors.BoundaryNorm(bounds, cmap.N)
ax = plt.subplot(111)

im = ax.scatter(x, y, c=t, s=S, marker='.',cmap='jet', edgecolors='none',norm=colors.LogNorm(vmin=1,vmax=1e5))
cbar = plt.colorbar(im)
cbar.ax.tick_params(labelsize=18)
plt.xlabel(r'$\Psi$',fontsize=24)
plt.ylabel(r'$I$',fontsize=24,rotation=0)
plt.axis([-0.5, 0.5, 0.33, 1.0])
plt.xticks([-0.5,-0.3,0.0,0.3, 0.5],[r'$-0.5$',r'$-0.3$',r'$0.0$',r'$0.3$',r'$0.5$'],fontsize=20)
plt.yticks([0.4,0.6,0.8,1.0],[r'$0.4$',r'$0.6$',r'$0.8$',r'$1.0$'],fontsize=20)
plt.savefig("Basins.jpg", dpi=500)
plt.show()
