import math

from scipy import integrate

import const

def integphi(mmq, mmuq, pphi, aaphi, T):
    a = 0.0
    L = const.L
    return integrate.quad(fphi, a, L, args=(mmq, mmuq, pphi, aaphi, T))[0]


def fphi(x, mq, muq, phi, aphi, T):
    E = math.sqrt(x * x + mq * mq)

    fer1 = math.exp(-(E - muq) / T)
    fer2 = math.exp(-(E + muq) / T)

    u1 = 1.0 + 3.0 * (phi + aphi * fer1) * fer1 + (fer1) ** 3
    u2 = 1.0 + 3.0 * (aphi + phi * fer2) * fer2 + (fer2) ** 3

    return x ** 2 * (fer1 / u1 + fer2 ** 2 / u2)
