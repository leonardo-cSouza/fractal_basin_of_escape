import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.gridspec as gridspec

plt.rc('text', usetex = True)
plt.rc('font', **{'family' : "sans-serif"})
params = {'text.latex.preamble' : [r'\usepackage{siunitx}', r'\usepackage{amsmath}']}
plt.rcParams.update(params)
plt.rcParams['axes.formatter.use_locale'] = True

# Import the file with the data
DATA = np.loadtxt("esta.dat") 

x = DATA[:,0]; y = DATA[:,1]; z = DATA[:,2]
S = 1.0

fig = plt.figure(figsize=(8, 6))
cmap = colors.ListedColormap(['blue', 'red'])
bounds=[1.0,2.0]
norm = colors.BoundaryNorm(bounds, cmap.N)

ax = plt.subplot(111)

im= ax.scatter(x, y, c=z, s=1.0, marker='.', edgecolors='none', cmap=cmap, norm=norm)
plt.xlabel(r'$\Psi/(2\pi)$',fontsize=24)
plt.ylabel(r'$I$',fontsize=24,rotation=0)
plt.axis([0.2, 0.4, 0.85, 0.95])
plt.xticks([0.2, 0.3, 0.4],[r'$0.2$',r'$0.3$',r'$0.4$'],fontsize=20)
plt.yticks(fontsize=20)
plt.savefig("unstable_col.jpg", dpi=300)
plt.show()