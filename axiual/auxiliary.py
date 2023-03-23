import math
from axiual import global_params_dict


def fermi1(en, mu, phi, aphi):
    fer1 = math.exp(-(en - mu) / global_params_dict['T'])
    u1 = 1.0 + 3.0 * (phi + aphi * fer1) * fer1 + fer1 ** 3

    return ((phi + 2.0 * aphi * fer1) * fer1 + fer1 ** 3) / u1


def fermi2(en, mu, phi, aphi, T):
    fer2 = math.exp(-(en + mu) / T)
    u2 = 1.0 + 3.0 * (aphi + phi * fer2) * fer2 + fer2 ** 3
    return ((aphi + 2.0 * phi * fer2) * fer2 + fer2 ** 3) / u2


def HeavisideTheta(arg):
    if (arg < 0.0):
        return 0.0
    else:
        return 1.0
