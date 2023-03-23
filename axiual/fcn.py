import numpy as np
from axiual import global_params_dict
from axiual.const import Gdiq, pi
from axiual.integral_A import integral_A
from axiual.integral_ImB import integral_ImB
from axiual.integral_ReB import integral_ReB


def fcn(x):
    f = np.zeros(2)
    ia_u = integral_A(global_params_dict['mq_u'], global_params_dict['mu_u'])
    # ia_d = integral_A(global_params_dict['mq_d'], global_params_dict['mu_d'])
    # ia_s = integral_A(global_params_dict['mq_s'], global_params_dict['mu_s'])
    # ia_s_min = integral_A(global_params_dict['mq_s'], - global_params_dict['mu_s'])
    print(global_params_dict['mq_u'], global_params_dict['phi'], global_params_dict['aphi'])
    ia_d_min = integral_A(global_params_dict['mq_d'], -global_params_dict['mu_d'])

    Gdiq_f = Gdiq / 4.0

    B0_ud = integral_ReB(x[0], global_params_dict['mq_u'], global_params_dict['mu_u'], global_params_dict['mq_d'],
                         -global_params_dict['mu_d'])

    # B0_us = integral_ReB(x[0], global_params_dict['mq_u'], global_params_dict['mu_u'], global_params_dict['mq_s'], -global_params_dict['mu_s'])

    Im_B0ud = integral_ImB(x[0], global_params_dict['mq_u'], global_params_dict['mu_u'], global_params_dict['mq_d'],
                           -global_params_dict['mu_d'])
    # Im_B0_us = integral_ImB(x[0], global_params_dict['mq_u'], global_params_dict['mu_u'], global_params_dict['mq_s'],
    #                        -global_params_dict['mu_s'])

    modul_du = B0_ud ** 2 + Im_B0ud ** 2

    #
    # Pol_ds = - 2.0 / pi ** 2 * (ia_u + ia_s_min
    #                             + (mq_u ** 2 + mq_s ** 2 - 4.0 * mq_u * mq_s
    #                                - (x[0] + mu_u + mu_s) ** 2) * B0_us)

    f[0] = Gdiq_f * (x[0] + global_params_dict['mu_u'] + global_params_dict['mu_d']) * x[1]
    + (1 - Gdiq_f * (ia_u + ia_d_min)) * Im_B0ud / modul_du

    f[1] = Gdiq_f * ((global_params_dict['mu_u'] ** 2 + global_params_dict['mu_d'] ** 2 - 4.0 * global_params_dict[
        'mu_u'] * global_params_dict['mu_d']
                      - (x[0] + global_params_dict['mu_u'] + global_params_dict['mu_d']) ** 2) + 1.0 / 4.0 * x[1] ** 2)
    + (1 - Gdiq_f * (ia_u + ia_d_min)) * B0_ud / modul_du

    return f
