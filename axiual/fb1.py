import math

from auxiliary import fermi1
from axiual import global_params_dict


def f1_1(x):
    f1 = fermi1(x, global_params_dict['mu_q1'], global_params_dict['phi'], global_params_dict['aphi'])
    p = x ** 2 - global_params_dict['mq1'] ** 2
    if p < 0.0:
        p = 0.0
    return math.sqrt(p) * f1 / 2.0 / global_params_dict['lam'] / (x - global_params_dict['E01'])


def f1_2(x):
    f1 = fermi1(x, global_params_dict['mu_q1'], global_params_dict['phi'], global_params_dict['aphi'])
    p = x ** 2 - global_params_dict['mq1']  ** 2
    if p < 0.0:
        p = 0.0
    return math.sqrt(p) * f1 / 2.0 / global_params_dict['lam']
