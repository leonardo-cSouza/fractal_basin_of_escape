import numpy as np
from numba import vectorize, njit
import matplotlib.pyplot as plt
import matplotlib as mpl
import random
from joblib import Parallel, delayed

q1, q2, q3 = 5, -6.3, 6.3
e1, e2, e3 = 10.7, -15.8, 4.13
v1, v2, v3, v4 = -9.867, 17.47, 10.1, -9
M, L, alpha, gamma = 15, 6, 1.83e-2, -9.16e-1
Niter = int(5e3)
beta = 0.04
parameters = np.array([q1, q2, q3, e1, e2, e3, v1, v2, v3, v4, M, L, alpha, beta, gamma])

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
    plt.rc('font', **{'family' : "sans-serif"})
    plt.rcParams["mathtext.fontset"] = 'dejavuserif'
    cmap = 'nipy_spectral'
    mpl.rcParams['axes.linewidth'] = 1.3 #set the value globally

    return cmap

@njit
def dig(psi0, I0, parameters, N):
    q1, q2, q3, e1, e2, e3, v1, v2, v3, v4, M, L, alpha, kappa, gamma = parameters
    
    pi = np.pi
    psi = psi0
    I = I0
    S = 0.0
    WB0 = 0.0
    WB1 = 0.0 
    for ii in range(1, N):
        u = float(ii)/float(N)
        S += np.exp(-((u*(1 - u))**(-1)))
    for i in range(N+1):
        u = float(i)/float(N)
        if (u > 0 and u < 1):
            w = np.exp(-((u*(1 - u))**(-1)))/S
            WB0 = WB0 + w*np.cos(psi)
            # Iterates the map
        I = I + kappa*np.sin(2*pi*psi)
        q_factor = q1 + q2*I**2 + q3*I**3
        v_parallel = v1 + v2*np.tanh(v3*I + v4)
        E_r = e1*I + e2*np.sqrt(abs(I)) + e3
        psi = (psi + alpha*v_parallel*(M/q_factor - L) + gamma*E_r/np.sqrt(abs(I)) + 0.5) % 1.0 - 0.5
    for i in range(N+1):
        u = float(i)/float(N)
        if(u > 0 and u < 1):
            w = np.exp(-((u*(1 - u))**(-1)))/S
            WB1 = WB1 + w*np.cos(psi)
            # Iterates the map           
        I = I + kappa*np.sin(2*pi*psi)
        q_factor = q1 + q2*I**2 + q3*I**3
        v_parallel = v1 + v2*np.tanh(v3*I + v4)
        E_r = e1*I + e2*np.sqrt(abs(I)) + e3
        psi = (psi + alpha*v_parallel*(M/q_factor - L) + gamma*E_r/np.sqrt(abs(I)) + 0.5) % 1.0 - 0.5

    d = -np.log10(np.abs(WB0 - WB1))
    
    return d

# Plot the probability distribution of the finite-time Lyapunov exponent

L = 1000
N = int(2e6)
psi = np.linspace(-0.5, 0.5, L, endpoint=True)
I = np.linspace(0.5, 1, L, endpoint=True)
psi, I = np.meshgrid(psi, I)

digs = np.zeros(L+1)
digs = Parallel(n_jobs=-1)(delayed(dig)(psi[i, j], I[i, j], parameters, N) for i in range(L) for j in range(L))

plot_params()
fig, ax = plt.subplots()
hm = ax.pcolor(psi, I, np.array(digs).reshape((L, L)), vmin=0, vmax=16, cmap="nipy_spectral", shading="auto")
plt.ylabel('$I$')
plt.xlabel('$\Psi$')
plt.tight_layout()
#plt.savefig("grid_dig_liul.jpg", dpi=300)
plt.show()
