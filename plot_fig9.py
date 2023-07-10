import numpy as np
from pylab import *
from matplotlib import colors

plt.rc('text', usetex = True)
plt.rc('font', **{'family' : "sans-serif"})
params = {'text.latex.preamble' : [r'\usepackage{siunitx}', r'\usepackage{amsmath}']}
plt.rcParams.update(params)
plt.rcParams['axes.formatter.use_locale'] = True

# import the data
DATA = np.loadtxt("wada_data.dat") 

x = DATA[:,0] 
y = DATA[:,1] 
z = DATA[:,2] 
t = DATA[:,3] 

xmin = x.min()
xmax = x.max()
ymin = y.min()
ymax = y.max()

S = 1.0

fig = plt.figure(figsize=(8, 6))
cmap = colors.ListedColormap(['lime', 'red', 'blue'])
bounds=[0, np.pi/2.0, np.pi,2.0*np.pi]
norm = colors.BoundaryNorm(bounds, cmap.N)
ax = plt.subplot(111)

im = ax.scatter(x, y, c=z, s=S, marker='.', edgecolors='none', cmap=cmap, norm=norm)
var = np.loadtxt('instavel_1.dat')
'''
#add this to plot the manifold with the basins 
x2 = var[:,0]; y2 = var[:,1]
for jj in range(len(y2)):
    if (y2[jj]>0.89 and y2[jj]<0.9):
        plt.scatter(x2[jj],y2[jj], color='yellow', s=1.0e-1)
'''
        plt.axis([-0.5, 0.5, 0.3, 1.005])
plt.xlabel(r'$\Psi/(2\pi)$',fontsize=24)
plt.ylabel(r'$I$',fontsize=24,rotation=0)

plt.xticks([-0.5,-0.3,0.0,0.3, 0.5],[r'$-0.5$',r'$-0.3$',r'$0.0$',r'$0.3$',r'$0.5$'],fontsize=20)
plt.yticks([0.4,0.6,0.8,1.0],[r'$0.4$',r'$0.6$',r'$0.8$',r'$1.0$'],fontsize=20)

# Make the square

plt.plot([-0.2,-0.2], [0.85,0.92], linestyle= '-', color = 'black',lw = 1)
plt.plot([0.0,0.0], [0.85,0.92], linestyle= '-', color = 'black',lw = 1)
plt.plot([-0.2,0.0], [0.85,0.85], linestyle= '-', color = 'black',lw = 1)
plt.plot([-0.2,0.0], [0.92,0.92], linestyle= '-', color = 'black',lw = 1)

plt.plot([-0.5,-1.0/6.0], [1.0,1.0], linestyle= '-', color = 'green',lw = 2.5)
plt.plot([-1.0/6.0,1.0/6.0], [1.0,1.0], linestyle= '-', color = 'blue',lw = 2.5)
plt.plot([1.0/6.0,0.5], [1.0,1.0], linestyle= '-', color = 'red',lw = 2.5)


plt.savefig("bacia-wada.jpg", dpi=300)

plt.show()
