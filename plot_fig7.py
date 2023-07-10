import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
import matplotlib.gridspec as gridspec
from matplotlib import rc
import matplotlib.cm as cm
import matplotlib as mpl
import locale

def plot_params(fontsize=14, tick_labelsize=17, axes_labelsize=20, legend_fontsize=14):
    """
    Update the parameters of the plot.

    Returns
    -------
    cmap : string
        The color map used in the colored plots.
    """
    plt.clf()
    plt.rc('font', size=fontsize)
    plt.rc('xtick', labelsize=tick_labelsize)
    plt.rc('ytick', labelsize=tick_labelsize)
    plt.rc('axes', labelsize=axes_labelsize)
    plt.rc('legend', fontsize=legend_fontsize)
    plt.rc('text', usetex = True)
    plt.rc('font', **{'family' : "sans-serif"})
    #plt.rc('font', **font)
    plt.rcParams["mathtext.fontset"] = "dejavuserif"
    cmap = 'nipy_spectral'
    mpl.rcParams['axes.linewidth'] = 1.3 #set the value globally

    return cmap

plt.rc('text', usetex = True)
plt.rc('font', **{'family' : "latin-modern-math"})
params = {'text.latex.preamble' : [r'\usepackage{siunitx}', r'\usepackage{amsmath}']}
plt.rcParams.update(params)
plt.rcParams['axes.formatter.use_locale'] = True

esta = np.loadtxt("stable.dat") 
inst = np.loadtxt("unstable.dat")
sela = np.loadtxt("saddle.dat")

x1 = esta[:,0]; y1 = esta[:,1] 

x2 = inst[:,0]; y2 = inst[:,1] 

x3 = sela[:,0]; y3 = sela[:,1] 

fig1 = plt.figure(figsize=(8, 7))
plot_params()
plt.scatter(x1, y1, color = 'red', s=0.5, edgecolor='none')
plt.xlabel(r'$\Psi$',fontsize=24)
plt.ylabel(r'$I$',fontsize=24,rotation=0)
plt.axis([-0.5, 0.5, 0.8, 1.0])
plt.xticks([-0.5,-0.3,0.0,0.3, 0.5],[r'$-0.5$',r'$-0.3$',r'$0.0$',r'$0.3$',r'$0.5$'],fontsize=20)
plt.yticks([0.8,0.9,1.0],[r'$0.8$',r'$0.9$',r'$1.0$'],fontsize=20)
plt.tight_layout()
plt.savefig("stable_01.jpg", dpi=500)    

fig2 = plt.figure(figsize=(8, 7))
plot_params()
plt.scatter(x2, y2, color = 'blue', s=.5, edgecolor='none')
plt.axis([-0.5, 0.5, 0.8, 1.0])
plt.xlabel(r'$\Psi$',fontsize=24)
plt.ylabel(r'$I$',fontsize=24,rotation=0)
plt.xticks([-0.5,-0.3,0.0,0.3, 0.5],[r'$-0.5$',r'$-0.3$',r'$0.0$',r'$0.3$',r'$0.5$'],fontsize=20)
plt.yticks([0.8,0.9,1.0],[r'$0.8$',r'$0.9$',r'$1.0$'],fontsize=20)
plt.tight_layout()
plt.savefig("unstable_01.jpg", dpi=500)    

fig3 = plt.figure(figsize=(8, 7))
plot_params()
plt.scatter(x3, y3, color = 'black', s=0.5, edgecolor='none')
plt.axis([-0.5, 0.5, 0.8, 1.0])
plt.xlabel(r'$\Psi$',fontsize=24)
plt.ylabel(r'$I$',fontsize=24,rotation=0)
plt.xticks([-0.5,-0.3,0.0,0.3, 0.5],[r'$-0.5$',r'$-0.3$',r'$0.0$',r'$0.3$',r'$0.5$'],fontsize=20)
plt.yticks([0.8,0.9,1.0],[r'$0.8$',r'$0.9$',r'$1.0$'],fontsize=20)
plt.tight_layout()
plt.savefig("saddle_01.jpg", dpi=500)    

plt.show()