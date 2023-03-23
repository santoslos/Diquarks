from axiual.auxiliary import fermi1, fermi2
import math
from axiual import global_params_dict


def fB3_1(x):
    f1 = fermi2(x, global_params_dict['mu_q2'], global_params_dict['phi'], global_params_dict['aphi'])
    p = x ** 2 - global_params_dict['mq2'] ** 2
    if p < 0.0:
        p = 0.0
    return math.sqrt(p) * f1 / 2.0 / global_params_dict['lam'] / (x + global_params_dict['E02'])


def fB3_2(x):
    f1 = fermi1(x, global_params_dict['mu_q2'], global_params_dict['phi'], global_params_dict['aphi'])
    p = x ** 2 - global_params_dict['mq2'] ** 2
    if p < 0.0:
        p = 0.0
    return math.sqrt(p) * f1 / 2.0 / global_params_dict['lam']
