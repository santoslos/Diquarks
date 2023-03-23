import math

from scipy import integrate

from auxiliary import fermi1, fermi2
import const


def integral_A(mmq, mmuq, pphi, aaphi, T):
    a = 0.0
    L = const.L
    print(integrate.quad(fA, a, L, args=(mmq, mmuq, pphi, aaphi, T))[1])
    return -4.0 * integrate.quad(fA, a, L, args=(mmq, mmuq, pphi, aaphi, T))[0]


def fA(x, mq, muq, phi, aphi, T):
    En = math.sqrt(x ** 2 + mq ** 2)

    b = 1.0 / T

    f = fermi1(x, mq, muq, phi, aphi, T)
    af = fermi2(x, mq, muq, phi, aphi, T)

    return x ** 2 / En * (1.0 - f - af)

