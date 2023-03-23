import math

from axiual import global_params_dict, fb1
from axiual.const import L
from scipy import integrate



def int_ReB1(mmq1, mmu_q1, mmq2, mmu_q2):
    global_params_dict['mq1'] = mmq1
    global_params_dict['mu_q1'] = mmu_q1
    global_params_dict['mq2'] = mmq2
    global_params_dict['mu_q2'] = mmu_q2

    if global_params_dict['mu_q2'] ** 2 - global_params_dict['mq2'] ** 2 + global_params_dict['mq1'] ** 2 > 0.0:
        a = math.sqrt(global_params_dict['mu_q2'] ** 2 - global_params_dict['mq2'] ** 2 + global_params_dict['mq1'] ** 2)
    else:
        a = global_params_dict['mq1']

    b = math.sqrt(L ** 2 + global_params_dict['mq1'] ** 2)

    if global_params_dict['E01'] < 0.0:
        return 4.0 * integrate.quad(fb1.f1_1, a, b, epsabs=0.0001, epsrel=0.0001)[0]
    else:
        c = global_params_dict['E01']
        if c == a:
            return 4.0 * integrate.quad(fb1.f1_1, a, b, points=[a, b], epsabs=0.0001, epsrel=0.0001)[0]
        else:
            return 4.0 * integrate.quad(fb1.f1_2, a, b, points=[c], epsabs=0.0001, epsrel=0.0001)[0]
