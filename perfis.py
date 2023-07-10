import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex = True)
plt.rc('font', **{'family' : "sans-serif"})
params = {'text.latex.preamble' : [r'\usepackage{siunitx}', r'\usepackage{amsmath}']}
plt.rcParams.update(params)
plt.rcParams['axes.formatter.use_locale'] = True

alpha = 3.570; beta = -7.921; gamma = 4.132
x = np.linspace(0.2, 1, 1000)
v = -9.867 + 17.478*np.tanh(10.1*x - 9.0) 
def q(x):
    if x <= 1:
        return 5.0 - 6.3*x**2 + 6.3*x**3
    else:
        return 5.0*x
q_values = [q(h) for h in x]

E = 3.0*alpha*x + 2.0*beta*np.sqrt(np.abs(x)) + gamma
  
fig1 = plt.figure(figsize=(8, 6))
plt.plot(x, q_values, color='black')
plt.xlabel('$I$', fontsize=24)
plt.ylabel('$q$', fontsize=24, rotation=0)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlim(0.2, 1.0)
plt.savefig("q_profile.jpg", dpi=300)

fig2 = plt.figure(figsize=(8, 6))
plt.plot(x, E, color='black')
plt.xlabel('$I$', fontsize=24)
#plt.ylabel('$E_r$', fontsize=24, rotation=0)
plt.text(0.1, -1.3,'$E_r$',fontsize=24)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlim(0.2, 1.0)
plt.savefig("E_profile.jpg", dpi=300)

fig3 = plt.figure(figsize=(8, 6))
plt.plot(x, v, color='black')
plt.xlabel('$I$', fontsize=24)
#plt.ylabel('$v_{\parallel}$', fontsize=24, rotation=0)
plt.text(0.1, -9.,'$v_{\parallel}$',fontsize=24)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlim(0.2, 1.0)
plt.savefig("v_profile.jpg", dpi=300)

plt.show()