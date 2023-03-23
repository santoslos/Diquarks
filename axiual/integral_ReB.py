from axiual import global_params_dict
from axiual.int_ReB1 import int_ReB1
from axiual.int_ReB2 import int_ReB2
from axiual.int_ReB3 import int_ReB3
from axiual.int_ReB4 import int_ReB4


def integral_ReB(mmes, mmq1, mmu_q1, mmq2, mmu_q2):
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

    Re_B1 = int_ReB1(global_params_dict['mq1'], global_params_dict['mu_q1'], global_params_dict['mq2'],
                     global_params_dict['mu_q2'])

    Re_B2 = int_ReB2(global_params_dict['mq1'], global_params_dict['mu_q1'], global_params_dict['mq2'],
                     global_params_dict['mu_q2'])

    Re_B3 = int_ReB3(global_params_dict['mq1'], global_params_dict['mu_q1'], global_params_dict['mq2'],
                     global_params_dict['mu_q2'])

    Re_B4 = int_ReB4(global_params_dict['mq1'], global_params_dict['mu_q1'], global_params_dict['mq2'],
                     global_params_dict['mu_q2'])

    return -Re_B1 - Re_B2 + Re_B3 + Re_B4
