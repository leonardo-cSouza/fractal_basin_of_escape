from cProfile import label
from telnetlib import theNULL
import numpy as np
from pylab import *
from matplotlib import colors
import matplotlib.gridspec as gridspec


rc('text', usetex = True)
rc('font', **{'family' : "sans-serif"})
params = {'text.latex.preamble' : [r'\usepackage{siunitx}', r'\usepackage{amsmath}']}
plt.rcParams.update(params)
plt.rcParams['axes.formatter.use_locale'] = True

# importa os dados do arquivo
DATA = np.loadtxt("W2W3_mt.dat") 
histo = np.loadtxt("HISTOGRAM_mt.dat")
xhist = histo[:,0]
yhist = histo[:,1]

x = DATA[:,0] # primeira coluna
y = DATA[:,1] # segunda coluna
z = DATA[:,2] # terceira coluna
#t = DATA[:,3] # terceira coluna


cmap = colors.ListedColormap(['lime', 'red','blue'])
bounds=[0,1, 2, 3]
norm = colors.BoundaryNorm(bounds, cmap.N)

fig1 = plt.figure(figsize=(8, 6))
plt.scatter(x, y,label='$W_2$')
plt.scatter(x, z, label='$W_3$')
plt.legend(fontsize=16)
plt.yticks([0,0.25,0.5,0.75,1.0],[r'$0.0$',r'$0.25$',r'$0.5$',r'$0.75$',r'$1.0$'],fontsize=20)
plt.xticks([0,5, 10, 15,20], fontsize=20)
plt.xlabel(r'$q$',fontsize=24)

DATA = np.loadtxt("grid.dat") 

x = DATA[:,0] # primeira coluna
y = DATA[:,1] # segunda coluna'
z = DATA[:,2] # terceira coluna

S = 0.8

fig2 = plt.figure(figsize=(8, 6))

cmap = colors.ListedColormap(['black','red','green'])
bounds=[0, 1.2, 2.2, 3.2]
norm = colors.BoundaryNorm(bounds, cmap.N)

# Grafico 3D com scatter

ax = plt.subplot(111)
plt.axis([0, 1e3, 0, 1e3])
im = ax.scatter(x, y, c=z, s=0.7, marker='.', edgecolors='none', cmap=cmap, norm=norm)

plt.xticks([0,250,500,750,1000],[r'$-0.5$',r'$-0.25$',r'$0.0$',r'$0.25$',r'$0.5$'],fontsize=20)
plt.yticks([200, 480,720,1.0e3],[r'$0.4$',r'$0.6$',r'$0.8$',r'$1.0$'],fontsize=20)
plt.xlabel(r'$\Psi/(2\pi)$',fontsize=24)
plt.ylabel(r'$I$',fontsize=24,rotation=0)
#savefig("grid.jpg", dpi=300)

fig3 = plt.figure(figsize=(8, 6))
#plt.axis([0, 20, 0, 45.5])
plt.bar(xhist,yhist)
plt.xticks([0,5, 10, 15,20], fontsize=20)
plt.yscale('log')
plt.yticks(fontsize=20)
#plt.text(-0.1, 46.0,r'$\times 10^3$',fontsize=14)
plt.xlabel(r'$q$',fontsize=20)
plt.ylabel(r'$\log{N}$',fontsize=20)
savefig("hist.jpg", dpi=300)
show()
