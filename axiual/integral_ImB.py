import math

from axiual import global_params_dict
from axiual.auxiliary import HeavisideTheta, fermi1, fermi2
from axiual.int_ReB1 import int_ReB1
from axiual.int_ReB2 import int_ReB2
from axiual.int_ReB3 import int_ReB3
from axiual.int_ReB4 import int_ReB4
from const import L, pi


def integral_ImB(mmes, mmq1, mmu_q1, mmq2, mmu_q2):
    Im_B1, Im_B2, Im_B3, Im_B4 = 0.0, 0.0, 0.0, 0.0
    mes = mmes
    global_params_dict['mq1'] = mmq1
    global_params_dict['mu_q1'] = mmu_q1
    global_params_dict['mq2'] = mmq2
    global_params_dict['mu_q2'] = mmu_q2
    global_params_dict['lam'] = mes - global_params_dict['mu_q2'] + global_params_dict['mu_q1']
    global_params_dict['E01'] = (global_params_dict['lam'] ** 2 + global_params_dict['mq1'] ** 2 - global_params_dict[
        'mq2'] ** 2) / 2.0 / global_params_dict['lam']
    global_params_dict['E02'] = (global_params_dict['lam'] ** 2 - global_params_dict['mq1'] ** 2 + global_params_dict[
        'mq2'] ** 2) / 2.0 / global_params_dict['lam']

    Theta = HeavisideTheta(global_params_dict['E01'] - global_params_dict['mq1']) * HeavisideTheta(
        math.sqrt(L ** 2 + global_params_dict['mq1'] ** 2) - global_params_dict['E01'])

    if Theta > 0.0:
        Im_B1 = 4.0 * pi * math.sqrt(global_params_dict['E01'] ** 2 - global_params_dict['mq1'] ** 2) * fermi1(
            global_params_dict['E01'],
            global_params_dict['mu_q1'],
            global_params_dict['phi'],
            global_params_dict['aphi']) * Theta / 2.0 / global_params_dict['lam']

    Theta = HeavisideTheta(-global_params_dict['E01'] - global_params_dict['mq1']) * HeavisideTheta(
        math.sqrt(L ** 2 + global_params_dict['mq1'] ** 2) + global_params_dict['E01'])

    if Theta > 0.0:
        Im_B2 = 4.0 * pi * math.sqrt(global_params_dict['E01'] ** 2 - global_params_dict['mq1'] ** 2) * fermi1(
            global_params_dict['E01'],
            global_params_dict['mu_q1'],
            global_params_dict['phi'],
            global_params_dict['aphi']) * Theta / 2.0 / global_params_dict['lam']

    Theta = HeavisideTheta(-global_params_dict['E02'] - global_params_dict['mq2']) * HeavisideTheta(
        math.sqrt(L ** 2 + global_params_dict['mq2'] ** 2) + global_params_dict['E02'])

    if Theta > 0.0:
        Im_B3 = 4.0 * pi * math.sqrt(global_params_dict['E02'] ** 2 - global_params_dict['mq2'] ** 2) * (
                1 - fermi2(global_params_dict['E02'],
                           global_params_dict['mu_q2'],
                           global_params_dict['phi'],
                           global_params_dict['aphi'])) * Theta / 2.0 / global_params_dict['lam']

    Theta = HeavisideTheta(global_params_dict['E02'] - global_params_dict['mq2']) * HeavisideTheta(
        math.sqrt(L ** 2 + global_params_dict['mq2'] ** 2) - global_params_dict['E02'])

    if Theta > 0.0:
        Im_B4 = 4.0 * pi * math.sqrt(global_params_dict['E02'] ** 2 - global_params_dict['mq2'] ** 2) * (
                1 - fermi2(global_params_dict['E02'],
                           global_params_dict['mu_q2'],
                           global_params_dict['phi'],
                           global_params_dict['aphi'])) * Theta / 2.0 / global_params_dict['lam']

    return Im_B1 + Im_B2 - Im_B3 - Im_B4
