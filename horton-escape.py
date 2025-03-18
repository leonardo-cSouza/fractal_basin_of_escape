#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 12:08:39 2024

@author: Leonardo Costa de Souza
@mail: lcsouza@if.usp.br
"""

from numba import njit, prange
from multiprocessing import Pool
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
import random
import matplotlib as mpl
import numpy as np
from matplotlib.colors import ListedColormap
import matplotlib.patches as patches

# Plots parameters
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
    font = {'family' : 'stix'}
    plt.rc('font', **font)
    plt.rcParams["mathtext.fontset"] = "stix"
    cmap = 'nipy_spectral'
    mpl.rcParams['axes.linewidth'] = 1.3 #set the value globally

    return cmap

#%%
# Define the parameters of the map
q1 = 5.0; q2 = -6.3; q3 = 6.e3
e1 = 10.7; e2 = -15.8; e3 = 4.13
v1 = -9.867; v2 = 17.47; v3 = 10.1; v4 = -9.0 
M = 15; L = 6; mu = 1.83e-2; rho = -9.16e-1

@njit
def map(x: float, y: float, eps: float) -> tuple:    
    yn = y + eps*np.sin(2.0*np.pi*x)
    q = q1 + q2*yn**2 + q3*yn**3
    v_para= v1 + v2*np.tanh(v3*yn + v4)
    E = e1*yn + e2*np.sqrt(np.abs(yn)) + e3
    xn = np.mod(x + mu*v_para*(M/q - L) + rho*E/np.sqrt(np.abs(yn)) + 0.5, 1.0) - 0.5

    return xn, yn

#%%
@njit
def time_series(x0: np.float64, y0: np.float64, eps: np.float64, N: np.int32) -> np.ndarray:
    u = np.zeros((N + 1, 2))

    x = x0
    y = y0

    u[0, 0] = x0
    u[0, 1] = y0
    for j in range(1, N):
        
        x, y = map(x, y, eps)

        u[j, 0] = x
        u[j, 1] = y

    return u
#%%
##
p = 0.1
##
#%% PLot o espaço de fase


random.seed(14)
n_ic = 75
plot_params()
fig = plt.figure(figsize=(8, 6))
for i in range(n_ic):
    x = random.random()
    y = random.random()
    x = -0.5 + 1.0*x
    y = 0.4 + (1 - 0.4)*y
    time = time_series(x, y, p, int(5e4))
    plt.plot(time[:, 0], time[:, 1], "kx", markersize=.025)

plt.xlim(-0.5, 0.5)
plt.ylim(0.375, 1)
plt.xlabel(r"$\Psi$"), plt.ylabel(r"$I$");
plt.tight_layout()
plt.savefig("phase_space")

#%%
@njit
def escape_time_and_basin(x0: np.float64, y0: np.float64, p: np.float64, y_esc: np.float64, Nmax: np.int32) -> np.ndarray:
    x = x0
    y = y0
    basin = 0  # Initialize basin to 0 (no escape)
    escape_time = Nmax  # Initialize escape_time to Nmax (no escape)

    for i in range(1, Nmax + 1):
        x, y = map(x, y, p)
        if y > y_esc:  # Check if y has crossed the escape threshold
            escape_time = i
            if -0.5 <= x <= 0.0:
                basin = 1  # Basin 1 if x is in (0, pi)
            elif 0.0 < x <= 0.5:
                basin = 2  # Basin 2 if x is in (pi, 2*pi)
            break
    
    return np.array([basin, escape_time])
#%%

# Grid setup
L = 2000
x0 = np.linspace(-0.5, 0.5, L, endpoint=True)#np.linspace(0, 2.0 * np.pi, L, endpoint=True)
y0 = np.linspace(0.375, 1.0, L, endpoint=True)#np.linspace(0.0, 0.75, L, endpoint=True)
x0, y0 = np.meshgrid(x0, y0)

xz = np.linspace(-0.2, 0.0, L, endpoint=True)#np.linspace(0, 2.0 * np.pi, L, endpoint=True)
yz = np.linspace(0.8, 0.85, L, endpoint=True)#np.linspace(0.0, 0.75, L, endpoint=True)
xz, yz = np.meshgrid(xz, yz)

# Parameters
#p = 0.39
y_esc = 1.0
Nmax = int(2e3)
#%%
@njit(parallel=True)
def compute_basins_parallel(x0, y0, p, y_esc, Nmax):
    L = x0.shape[0]
    basins = np.zeros_like(x0, dtype=np.int32)
    time = np.zeros_like(x0, dtype=np.int32)
    for i in prange(L):
        for j in prange(L):
            basins[i, j] = escape_time_and_basin(x0[i, j], y0[i, j], p, y_esc, Nmax)[0]
            time[i, j] = escape_time_and_basin(x0[i, j], y0[i, j], p, y_esc, Nmax)[1]
    return basins, time
#%%
# Compute basins in parallel
eps = 0.1
basins, tempo = compute_basins_parallel(x0, y0, eps, y_esc, Nmax)
basins_z, tempo_z = compute_basins_parallel(xz, yz, eps, y_esc, Nmax)
#%%
# Plotting 
plot_params()
fig2 = plt.figure(figsize=(8, 6))

cmap = ListedColormap(['white', '#1c698a', '#ff904a'])


# Plot the grid, scatter based on the basin colors
plt.scatter(x0, y0, c=basins, cmap=cmap, marker='x', s=0.0005)  # Adjust marker size as needed
rect = patches.Rectangle((-0.2, 0.8), 0.2, 0.05, linewidth=2, edgecolor='black', facecolor='none')
plt.gca().add_patch(rect)
# Show colorbar for reference
#plt.colorbar()

plt.xlabel(r"$\Psi$"); plt.ylabel("$I$")
plt.xlim(-0.5, 0.5)
plt.ylim(0.36, 1)
plt.tight_layout()
path = "/home/leonardo/Documentos/PD/Figures/"
figname = path + "Escape-horton-eps=%.3f.png" % (eps)
plt.savefig(figname, dpi=300, format=figname[-3:], bbox_inches='tight', pad_inches=0.05)
print("Figure save in ", path)

#%%

# Plot the grid, scatter based on the basin colors
plot_params()
fig3 = plt.figure(figsize=(8, 6))
plt.scatter(xz, yz, c=basins_z, cmap=cmap, marker='x', s=0.0005)  # Adjust marker size as needed
plt.xlabel(r"$\Psi$"); plt.ylabel("$I$")
plt.xlim(-0.2, 0.0)
plt.ylim(0.8, 0.85)
plt.xticks([-0.2, -0.15, -0.1, -0.05, 0.0])
plt.tight_layout()
path = "/home/leonardo/Documentos/PD/Figures/"
figname = path + "Escape-horton-zoom-eps=%.3f.png" % (eps)
plt.savefig(figname, dpi=300, format=figname[-3:], bbox_inches='tight', pad_inches=0.05)
print("Figure save in ", path)

#%%
L = 1000
x0 = np.linspace(-0.5, 0.5, L, endpoint=True)#np.linspace(0, 2.0 * np.pi, L, endpoint=True)
y0 = np.linspace(0.7, 1.0, L, endpoint=True)#np.linspace(0.0, 0.75, L, endpoint=True)
x0, y0 = np.meshgrid(x0, y0)

# Compute basins in parallel
eps = 0.1

basins, tempo = compute_basins_parallel(x0, y0, eps, y_esc, Nmax)
# Parameters
itera_short = 6  # Number of short iterations

# Initialize lists to store the data
stable_data = []
unstable_data = []
saddle_data = []

@njit
def short_iteration(x0, y0, p, itera_short):
    x, y = x0, y0
    xsela, ysela = 0.0, 0.0

    for n in range(itera_short):
        x, y = map(x, y, p)

        if n == itera_short // 2:
            xsela, ysela = x, y

    return x, y, xsela, ysela

# Process points based on basins
for i in range(L):
    for j in range(L):
        if basins[i, j] != 0:  # Consider only points with non-zero basin value
            xi, yi = x0[i, j], y0[i, j]
            xn, yn, xsela, ysela = short_iteration(xi, yi, p, itera_short)

            stable = None
            unstable = None
            saddle = None

            if 0.0 < yn < 1.0 and 0.0 < yi < 1.0 and 0.0 < ysela < 1.0:
                stable = (xi, yi, 2.0)
                unstable = (xn, yn, 2.0)
                saddle = (xsela, ysela)

                if stable:
                    stable_data.append(stable)
                if unstable:
                    unstable_data.append(unstable)
                if saddle:
                    saddle_data.append(saddle)

# Save the results to files
np.savetxt('stable2.dat', stable_data, fmt='%.6f %.6f %.6f')
np.savetxt('unstable2.dat', unstable_data, fmt='%.6f %.6f %.6f')
np.savetxt('saddle2.dat', saddle_data, fmt='%.6f %.6f')

#%%
stable_manifold = np.array(stable_data)
unstable_manifold = np.array(unstable_data)
saddle_manifold = np.array(saddle_data)

#%%
plot_params()
fig, ax = plt.subplots(1, 3, figsize=(18, 8))
ax[0].scatter(stable_manifold[:, 0], stable_manifold[:, 1], color='RED', label='Stable', s=.005)
ax[1].scatter(unstable_manifold[:, 0], unstable_manifold[:, 1], color='navy', label='Stable', s=.005)
ax[2].scatter(saddle_manifold[:, 0], saddle_manifold[:, 1], color='black', label='Stable', s=.005)
ax[0].set_ylabel("$y$")
ax[0].set_xlabel("$x$")
ax[1].set_xlabel("$x$")
ax[2].set_xlabel("$x$")

ax[0].set_ylim(0.6, 1.0)
ax[1].set_ylim(0.6, 1.0)
ax[2].set_ylim(0.6, 1.0)
ax[0].set_xlim(-0.5, 0.5)
ax[1].set_xlim(-0.5, 0.5)
ax[2].set_xlim(-0.5, 0.5)
plt.tight_layout()
plt.savefig("manifolds.png", dpi=300)    
#%%
plot_params()
fig4 = plt.figure(figsize=(8, 6))
plt.scatter(stable_manifold[:, 0], stable_manifold[:, 1], color='RED', label='Stable', s=.005)
plt.xlabel(r"$\Psi$")
plt.ylabel(r"$I$")
plt.tight_layout()
plt.ylim(0.7, 1.0)
plt.xlim(-0.5, 0.5)
#%%
plot_params()
fig4 = plt.figure(figsize=(8, 6))
plt.scatter(unstable_manifold[:, 0], unstable_manifold[:, 1], color='navy', label='Stable', s=.01)
plt.xlabel(r"$\Psi$")
plt.ylabel(r"$I$")
plt.tight_layout()
plt.ylim(0.7, 1.0)
plt.xlim(-0.5, 0.5)

#%%
@njit
def dig(x0, y0, p, N):
    u = np.arange(1, N) / N
    S = np.sum(np.exp(-1 / (u * (1 - u))))

    ts = time_series(x0, y0, p, N)
    w = np.exp(-1 / (u * (1 - u))) / S
    WB0 = np.sum(w * np.cos(ts[1:-1, 0]))

    x0 = ts[-1, 0]
    y0 = ts[-1, 1]
    ts = time_series(x0, y0, p, N)
    w = np.exp(-1 / (u * (1 - u))) / S
    WB1 = np.sum(w * np.cos(ts[1:-1, 0]))

    return -np.log10(abs(WB0 - WB1))
#%%
@njit(parallel=True)
def compute_dig(x0, y0, p, Nmax):
    L = x0.shape[0]
    digs = np.zeros_like(x0, dtype=np.float64)

    for i in prange(L):
        for j in prange(L):
            digs[i, j] = dig(x0[i, j], y0[i, j], p, Nmax)
    
    return digs
#%%
# Example usage (assuming you have set up x0, y0, and other parameters)
# Compute digs in parallel
digs = compute_dig(x0, y0, 0.4, 100)
#%%
NN= int(1e4)

L = 100
x0 = np.linspace(0, 2.0 * np.pi, L, endpoint=True)
y0 = np.linspace(0.0, .2, L, endpoint=True)
x0, y0 = np.meshgrid(x0, y0)

# Parameters
param = 0.1
digs = compute_dig(x0, y0, param, NN)
#%%
# Create the heatmap
plt.figure(figsize=(10, 8))
#plt.imshow(digs, cmap='viridis', origin='lower', extent=(0, 2*np.pi, 0, 1))
plt.pcolor(x0, y0, digs.reshape(L, L), cmap='plasma',vmin = 0.0)
# Add color bar to indicate the scale of 'dig' values
plt.colorbar(label='dig')
# Add labels and title
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.show()

