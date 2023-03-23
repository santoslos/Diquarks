from scipy import integrate

from auxiliary import fermi1, fermi2
import const


def intmu(mq, muq, phi, aphi, T):
    a = 0.0
    L = const.L
    return integrate.quad(funcmu, a, L, args=(mq, muq, phi, aphi, T))[0]


def funcmu(x, mq, muq, phi, aphi, T):
    f1 = fermi1(x, mq, muq, phi, aphi, T)
    f2 = fermi2(x, mq, muq, phi, aphi, T)
    return x * x * (f1 - f2)
