import math

from axiual import global_params_dict
from axiual.const import L
from scipy import integrate
import fb4


def int_ReB4(mmq1, mmu_q1, mmq2, mmu_q2):
    global_params_dict['mq1'] = mmq1
    global_params_dict['mu_q1'] = mmu_q1
    global_params_dict['mq2'] = mmq2
    global_params_dict['mu_q2'] = mmu_q2
    a = global_params_dict['mu_q2']
    b = math.sqrt(L ** 2 + global_params_dict['mu_q2'] ** 2)

    if global_params_dict['E02'] < 0.0:
        return 4.0 * integrate.quad(fb4.fb_1, a, b, epsabs=0.0001, epsrel=0.0001)[0]
    else:
        c = global_params_dict['E02']
        if c == a:
            return 4.0 * integrate.quad(fb4.fb_1, a, b, points=[a, b], epsabs=0.0001, epsrel=0.0001)[0]
        else:
            return 4.0 * integrate.quad(fb4.fb_2, a, b, points=[c], epsabs=0.0001, epsrel=0.0001)[0]
