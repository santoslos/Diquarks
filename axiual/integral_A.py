import math

from scipy import integrate

from auxiliary import fermi2, fermi1
from axiual import global_params_dict
from axiual.const import L


def fmq(x):
    en = math.sqrt(global_params_dict['mq'] * global_params_dict['mq'] + x * x)
    return x ** 2 / en * (1.0 - fermi1(en, global_params_dict['mu_q'], global_params_dict['phi'], global_params_dict['aphi']) - fermi2(en,
                                                                                                               global_params_dict['mu_q'],
                                                                                                               global_params_dict['phi'],
                                                                                                               global_params_dict['aphi']))


def integral_A(mmq, mmu_q):
    global_params_dict['mq'] = mmq
    global_params_dict['mu_q'] = mmu_q
    print(mmq, mmu_q)
    if global_params_dict['mu_u'] ** 2 - global_params_dict['mq_u'] ** 2 > 0.0:
        a = math.sqrt(global_params_dict['mu_u'] ** 2 - global_params_dict['mq_u'] ** 2)
    else:
        a = 0.0

    return -4.0 * integrate.quad(fmq, a, L, epsabs=0.0001, epsrel=0.0001)[0]
