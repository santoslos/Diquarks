import numpy as np
import pandas as pd
import scipy

from axiual import global_params_dict
from axiual import fcn

x = np.zeros(2)
x[0] = 1.1
x[1] = 0.1

mueq_u_mass_phi = pd.read_csv('../mueq05Gv0GT_u_mass_phi_mu0.txt', header=None)
mueq_dmass = pd.read_csv('../mueq05Gv0GT_dmass_mu0.txt', header=None)
mueq_smass = pd.read_csv('../mueq05Gv0GT_smass_mu0.txt', header=None)

mq_u = mueq_u_mass_phi[2]
phi = mueq_u_mass_phi[3]
aphi = mueq_u_mass_phi[4]
mq_d = mueq_dmass[2]
mu_d = mueq_dmass[3]

T = mueq_smass[0]
mq_s = mueq_smass[2]
mu_s = mueq_smass[3]

for i in range(350):
    global_params_dict['mq_u'] = mq_u[i]
    global_params_dict['phi'] = phi[i]
    global_params_dict['aphi'] = aphi[i]
    global_params_dict['mq_d'] = mq_d[i]
    global_params_dict['mu_d'] = mu_d[i]
    global_params_dict['mq_s'] = mq_s[i]
    global_params_dict['mu_s'] = mu_s[i]
    global_params_dict['T'] = T[i]
    print(i)
    xguess = x.copy()
    x = scipy.optimize.root(fcn.fcn, xguess)['x']
    print(global_params_dict['T'], x[0])
