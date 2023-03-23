import math

from scipy import integrate

import const


def integaphi(mq, muq, phi, aphi, T):
    a = 0.0
    L = const.L
    return integrate.quad(faphi, a, L, args=(mq, muq, phi, aphi, T))[0]


def faphi(x, mq, muq, phi, aphi, T):
    E = math.sqrt(x * x + mq * mq)

    fer1 = math.exp(-(E - muq) / T)
    fer2 = math.exp(-(E + muq) / T)

    u1 = 1.0 + 3.0 * (phi + aphi * fer1) * fer1 + (fer1) ** 3
    u2 = 1.0 + 3.0 * (aphi + phi * fer2) * fer2 + (fer2) ** 3

    return x ** 2 * (fer1 ** 2 / u1 + fer2 / u2)
