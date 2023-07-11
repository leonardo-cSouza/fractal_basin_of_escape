import numpy as np
import matplotlib.pyplot as plt
import random

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

@njit
def map(psi0, I0, parameters, N):
    q1, q2, q3, e1, e2, e3, v1, v2, v3, v4, M, L, alpha, kappa, gamma = parameters
    
    u = np.zeros((N + 1, 2))

    pi = np.pi
    psi = psi0
    I = I0

    u[0, 0] = psi0
    u[0, 1] = I0

    for j in range(1, N):
        I = I + kappa*np.sin(2*pi*psi)
        q_factor = q1 + q2*I**2 + q3*I**3
        v_parallel = v1 + v2*np.tanh(v3*I + v4)
        E_r = e1*I + e2*np.sqrt(abs(I)) + e3
        psi = np.mod((psi + alpha*v_parallel*(M/q_factor - L) + gamma*E_r/np.sqrt(abs(I)) + 0.5), 1.0)- 0.5

        u[j, 0] = psi
        u[j, 1] = I

    return u

"""@vectorize(['float64(float64, float64, float64[:], int64)'],
           target='parallel',
           nopython=True)"""
q1, q2, q3 = 5, -6.3, 6.3
e1, e2, e3 = 10.7, -15.8, 4.13
v1, v2, v3, v4 = -9.867, 17.47, 10.1, -9
M, L, alpha, gamma = 15, 6, 1.83e-2, -9.16e-1
N = 1000
beta = 0.04
parameters = np.array([q1, q2, q3, e1, e2, e3, v1, v2, v3, v4, M, L, alpha, beta, gamma])
plot_params()
random.seed(14)
n_ic = 75
for i in range(n_ic):
    x = random.random()
    y = random.random()
    x = -0.5 + x
    y = 0.4 + (1 - 0.4)*y
    time_series = map(x, y, parameters, N)
    plt.plot(time_series[:, 0], time_series[:, 1], "kx", markersize=0.1)
    
plt.xlim(-0.5, 0.5)
plt.ylim(0.5, 1)
plt.xlabel("$\\Psi$"), plt.ylabel("$I$");
#plt.tight_layout()

plt.show()