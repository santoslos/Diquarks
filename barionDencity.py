from scipy import integrate

import const
from auxiliary import fermi1, fermi2


def funcden(x, mq, muq, phi, aphi, T):
    f1 = fermi1(x, mq, muq, phi, aphi, T)
    f2 = fermi2(x, mq, muq, phi, aphi, T)
    return x * x * (f1 - f2)


def intden(mq, muq, phi, aphi, T):
    a = 0.0
    L = const.L
    return integrate.quad(funcden, a, 2.0 * L, args=(mq, muq, phi, aphi, T))[0]
