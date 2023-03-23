import math


def fermi1(pp, mmq, mmuq, pphi, aaphi, T):
    p = pp
    mq = mmq
    muq = mmuq
    phi = pphi
    aphi = aaphi
    en = math.sqrt(mq ** 2 + p ** 2)
    fer1 = math.exp(-(en - muq) / T)
    u1 = 1.0 + 3.0 * (phi + aphi * fer1) * fer1 + (fer1) ** 3
    fermi1 = ((phi + 2.0 * aphi * fer1) * fer1 + (fer1) ** 3) / u1
    return fermi1


def fermi2(pp, mmq, mmuq, pphi, aaphi, T):
    p = pp
    mq = mmq
    muq = mmuq
    phi = pphi
    aphi = aaphi
    en = math.sqrt(mq ** 2 + p ** 2)
    fer2 = math.exp(-(en + muq) / T)
    u2 = 1.0 + 3.0 * (aphi + phi * fer2) * fer2 + (fer2) ** 3
    fermi2 = ((aphi + 2.0 * phi * fer2) * fer2 + (fer2) ** 3) / u2
    return fermi2
