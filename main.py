# This is a sample Python script.
import numpy as np
# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.


# See PyCharm help at https://www.jetbrains.com/help/pycharm/
import scipy.integrate as integrate
import numpy as np

import scipy.optimize

from IntegA import integral_A
from barionDencity import intden
from const import Nc, pi
from fcn import fcn

x = np.zeros(8)

x[0] = 0.37
x[1] = 0.37
x[2] = 0.55
x[3] = 0.6
x[4] = 0.6
x[5] = 0.0
x[6] = 0.0
x[7] = 0.0

T = 0.002

mueq05Gv0GT_u_mass_phi_mu0 = open('mueq05Gv0GT_u_mass_phi_mu0.txt', 'w+')
mueq05Gv0GT_dmass_mu0 = open('mueq05Gv0GT_dmass_mu0.txt', 'w+')
mueq05Gv0GT_smass_mu0 = open('mueq05Gv0GT_smass_mu0.txt', 'w+')
eqmu05Gv0GT_condens_mu0 = open('eqmu05Gv0GT_condens_mu0.txt', 'w+')


for i in range(350):
    T = 0.002
    T = T + 0.001 * i
    xguess = x.copy()
    x = scipy.optimize.root(fcn, xguess, args=(T))['x']
    mu_u = x[5]
    mu_d = x[6]
    mu_s = x[7]

    ro_u = intden(x[0], mu_u, x[3], x[4], T)
    ro_d = intden(x[1], mu_d, x[3], x[4], T)
    ro_s = intden(x[2], mu_s, x[3], x[4], T)

    ro_tot = (ro_u + ro_d + ro_s) / 3.0

    uu = Nc * x[0] / pi ** 2 * integral_A(x[0], mu_u, x[3], x[4], T)
    ss = Nc * x[2] / pi ** 2 * integral_A(x[2], mu_s, x[3], x[4], T)

    print(f'{T:.{3}f},         {x[0]:.{5}f},         {x[3]:.{5}f},         {x[2]:.{5}f},          {x[5]:.{5}f}')
    eqmu05Gv0GT_condens_mu0.write(f'{T:.{3}f}, {uu:.{5}f}, {ss:.{5}f}, {ss/uu}\n')
    mueq05Gv0GT_u_mass_phi_mu0.write(f'{T:.{3}f}, {ro_tot:.{5}f}, {x[0]:.{5}f}, {x[3]:.{5}f}, {x[4]:.{5}f}\n')
    mueq05Gv0GT_dmass_mu0.write(f'{T:.{3}f}, {ro_tot:.{5}f}, {x[1]:.{5}f}, {mu_d:.{5}f}\n')
    mueq05Gv0GT_smass_mu0.write(f'{T:.{3}f}, {ro_tot:.{5}f}, {x[2]:.{5}f}, {mu_s:.{5}f}\n')