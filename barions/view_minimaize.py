import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

mu = ['00', '01', '03', '05', '07', '09', '11', '13', '15', '17', '19', '21', '23', '25', '26', '27', '28', '29',
      '30', '31', '32', '33', '34', '35']
plt.rcParams['text.latex.preamble'] = r"\usepackage{bm} \usepackage{amsmath} '\boldmath"
# '36', '37']
Tsc_list = []
fig, ax = plt.subplots()
# lambda_bar = np.loadtxt('../../../NJL/Goshe/Goshe/04axial_diq_ud00.txt')
# T_l = lambda_bar[:, 0]
# l_b = lambda_bar[:, 1]
# mq_u2_mq_s = lambda_bar[:, 2]
# plt.plot(T_l, l_b, label='line 1')
# plt.plot(T_l, mq_u2_mq_s, label='line 2')
# plt.show()


# for i in mu:
#     scalar = np.loadtxt(f'../../Gosha/Barions/Goshe/04axial_diq_ss{i}.txt')
#     # # if i in ['00', '01', '03', '05', '07',
#     #                                                                                   '09', '11'] else np.loadtxt(
#     #     f'../../Gosha/Barions/Goshe/04axial_diq{i}.txt')
#     m_u = np.loadtxt(f'../../Gosha/Barions/Goshe/04mueq_u_mass_phi_mu{i}.txt')
#     m_s = np.loadtxt(f'../../Gosha/Barions/Goshe/04mueq_smass_mu{i}.txt')
#     T_p = scalar[:, 0]
#     print(len(T_p))
#     m_ud_sc = scalar[0:280, 2] #if i in ['00', '01', '03', '05', '07', '09', '11'] else scalar[0:280, 2]
#     mq_u = 2 * m_u[0:280, 2]
#     mq_s = m_s[0:280, 2]
#     porog = mq_u + mq_s
#     Tsc_list.append(T_p[np.argmin(np.absolute(2*mq_s - porog))])
# # # print(Td_List)
# # print(T_lambda_bar_list)
#
# mu_2 = [0.00, 0.01, 0.03, 0.05, 0.07, 0.09, 0.11, 0.13, 0.15, 0.17, 0.19, 0.21, 0.23, 0.25, 0.26, 0.27, 0.28, 0.29,
#         0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37]
# # np.savetxt('T_mu_lambda.txt', np.c_[mu_2, T_lambda_bar_list], delimiter='    ')
# # np.savetxt('T_mu_protons.txt', np.c_[mu_2, Td_List], delimiter='    ')
# # np.savetxt('T_mu_psi.txt', np.c_[mu_2, T_psi_list], delimiter='    ')
#
# plt.scatter(mu, Tsc_list, c="g", label="axial_ss")
# # plt.scatter(mu_2, T_lambda_bar_list, c="r", label="lambda", marker="^")
# # plt.scatter(mu_2, T_psi_list, c="b", label="psi", marker="*")

#
axial_ss = np.loadtxt(f'../../Gosha/Barions/Goshe/04axial_diq_ss00.txt')
axial_us = np.loadtxt(f'../../Gosha/Barions/Goshe/04axial_diq_us00.txt')
axial_ud = np.loadtxt(f'../../Gosha/Barions/Goshe/04axial_diq_ud00.txt')

scalar_ud = np.loadtxt(f'../../Gosha/Barions/Goshe/04scalar_diq_ud00.txt')
scalar_us = np.loadtxt(f'../../Gosha/Barions/Goshe/04scalar_diq_us00.txt')

m_u = np.loadtxt(f'../../Gosha/Barions/Goshe/04mueq_u_mass_phi_mu00.txt')
m_s = np.loadtxt(f'../../Gosha/Barions/Goshe/04mueq_smass_mu00.txt')

T_axial_ss = axial_ss[:300, 0]
M_axial_ss = axial_ss[:300, 1]

T_axial_us = axial_us[:300, 0]
M_axial_us = axial_us[:300, 1]

T_axial_ud = axial_ud[:300, 0]
M_axial_ud = axial_ud[:300, 1]

T_scalar_ud = scalar_ud[:300, 0]
M_scalar_ud = scalar_ud[:300, 1]

T_scalar_us = scalar_us[:300, 0]
M_scalar_us = scalar_us[:300, 1]

T_mu = m_u[:300, 0]
mq_u = m_u[:300, 2]

mq_s = m_s[:300, 2]

ax.plot(T_axial_ss, M_axial_ss, label='[ss]', linewidth=3)
ax.plot(T_axial_us, M_axial_us, label='[us]', linewidth=3, color='blue')
ax.plot(T_axial_ud, M_axial_ud, label='[ud]', linewidth=3, color='brown')
ax.plot(T_scalar_ud, M_scalar_ud, label='(ud)', linewidth=3, color='black')
ax.plot(T_scalar_us, M_scalar_us, label='(us)', linewidth=3, color='red')

ax.plot(T_mu, 2 * mq_u, '--', label='$2m_u', linewidth=2, color='black')
ax.plot(T_mu, mq_s + mq_u, '-.', label='$m_u + $m_s', linewidth=2, color='black')
ax.plot(T_mu, 2 * mq_s, ':', label='2$m_s', linewidth=2)
plt.text(0.03, 1.04, '[ss]', fontsize=18)
plt.text(0.03, 0.92, '[us]', fontsize=18)
plt.text(0.273, 0.83, '[ud]', fontsize=18)
plt.text(0.03, 0.60, '(ud)', fontsize=18)
plt.text(0.25, 0.64, '(us)', fontsize=18)

plt.text(0.23, 0.02, r'${2m_u}$', fontsize=18)
plt.text(0.235, 0.3, r'${m_u + m_s}$', fontsize=18)
plt.text(0.153, 1.1, r'${2m_s}$', fontsize=18)
ax.set(xlim=[0.003, 0.3], ylim=[0, 1.2])
ax.tick_params(which='major', direction='in', width=2, length=6)
for axis in ['top', 'bottom', 'left', 'right']:
    ax.spines[axis].set_linewidth(2)
plt.xlabel("Т, GeV", fontsize=18)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.ylabel("masses, GeV", fontsize=18)
plt.show()
# fig.savefig('Diquarks0.eps')
#
# proton = np.loadtxt(f'../../Gosha/Barions/Goshe/04proton00.txt')
# lambda_bar = np.loadtxt(f'../../Gosha/Barions/Goshe/04lambda_bar00.txt')
# psi = np.loadtxt(f'../../Gosha/Barions/Goshe/04psi_bar00.txt')
#
# T_p = proton[:, 0]
# M_p = proton[:182, 1]
# mu3 = proton[:, 2]
#
# T_l = lambda_bar[:, 0]
# l_b = lambda_bar[:197, 1]
# mq_u2_mq_s = lambda_bar[:, 2]
#
# T_psi = psi[:, 0]
# psi_bar = psi[:239, 1]
# mq_u2_mq_s_psi = psi[:, 2]
#
# ax.plot(T_p[:182], M_p, label='$p$', linewidth=3, color='blue')
# ax.plot(T_l[:197], l_b, label='$\Lambda^0$', linewidth=3, color='brown')
# ax.plot(T_psi[:239], psi_bar, label='$\Xi^-$', linewidth=3, color='green')
#
# ax.plot(T_p, mu3, '--', label='$2m_u$', linewidth=2, color='blue')
# ax.plot(T_psi, mq_u2_mq_s_psi, '-.', label='$2m_u$ + $m_s$', linewidth=2, color='green')
# ax.plot(T_l, mq_u2_mq_s, ':', label='$m_u$ + $2m_s$', linewidth=2, color='brown')
#
# plt.text(0.21, 1.2, '$m_u + 2m_s$', fontsize=18)
# plt.text(0.23, 0.3, '$2m_u + m_s$', fontsize=18)
# plt.text(0.21, 0.11, '$2m_u$ ', fontsize=18)
# plt.text(0.19, 0.9, '$\Lambda^0$', fontsize=20)
# plt.text(0.03, 1.29, '$\Xi^-$ ', fontsize=20,)
# plt.text(0.03, 0.92, '$p$ ', fontsize=20)
# ax.set(xlim=[0.003, 0.3], ylim=[0, 1.55])
# ax.tick_params( which='major', direction='in', width=2, length=6)
# for axis in ['top','bottom','left','right']:
#     ax.spines[axis].set_linewidth(2)
# plt.xlabel("Т, GeV", fontsize=18 )
# plt.xticks(fontsize=15)
# plt.yticks(fontsize=15)
# plt.ylabel("masses, GeV", fontsize=18)
# plt.show()
# fig.savefig('Barions0.eps')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

mu = ['00', '01', '03', '05', '07', '09', '11', '13', '15', '17', '19', '21', '23', '25', '26', '27', '28', '29',
      '30', '31', '32', '33', '34', '35']
plt.rcParams['text.latex.preamble'] = r"\usepackage{bm} \usepackage{amsmath} '\boldmath"
# '36', '37']
Tsc_list = []
fig, ax = plt.subplots()
# lambda_bar = np.loadtxt('../../../NJL/Goshe/Goshe/04axial_diq_ud00.txt')
# T_l = lambda_bar[:, 0]
# l_b = lambda_bar[:, 1]
# mq_u2_mq_s = lambda_bar[:, 2]
# plt.plot(T_l, l_b, label='line 1')
# plt.plot(T_l, mq_u2_mq_s, label='line 2')
# plt.show()


# for i in mu:
#     scalar = np.loadtxt(f'../../Gosha/Barions/Goshe/04axial_diq_ss{i}.txt')
#     # # if i in ['00', '01', '03', '05', '07',
#     #                                                                                   '09', '11'] else np.loadtxt(
#     #     f'../../Gosha/Barions/Goshe/04axial_diq{i}.txt')
#     m_u = np.loadtxt(f'../../Gosha/Barions/Goshe/04mueq_u_mass_phi_mu{i}.txt')
#     m_s = np.loadtxt(f'../../Gosha/Barions/Goshe/04mueq_smass_mu{i}.txt')
#     T_p = scalar[:, 0]
#     print(len(T_p))
#     m_ud_sc = scalar[0:280, 2] #if i in ['00', '01', '03', '05', '07', '09', '11'] else scalar[0:280, 2]
#     mq_u = 2 * m_u[0:280, 2]
#     mq_s = m_s[0:280, 2]
#     porog = mq_u + mq_s
#     Tsc_list.append(T_p[np.argmin(np.absolute(2*mq_s - porog))])
# # # print(Td_List)
# # print(T_lambda_bar_list)
#
# mu_2 = [0.00, 0.01, 0.03, 0.05, 0.07, 0.09, 0.11, 0.13, 0.15, 0.17, 0.19, 0.21, 0.23, 0.25, 0.26, 0.27, 0.28, 0.29,
#         0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37]
# # np.savetxt('T_mu_lambda.txt', np.c_[mu_2, T_lambda_bar_list], delimiter='    ')
# # np.savetxt('T_mu_protons.txt', np.c_[mu_2, Td_List], delimiter='    ')
# # np.savetxt('T_mu_psi.txt', np.c_[mu_2, T_psi_list], delimiter='    ')
#
# plt.scatter(mu, Tsc_list, c="g", label="axial_ss")
# # plt.scatter(mu_2, T_lambda_bar_list, c="r", label="lambda", marker="^")
# # plt.scatter(mu_2, T_psi_list, c="b", label="psi", marker="*")

#
axial_ss = np.loadtxt(f'../../Gosha/Barions/Goshe/04axial_diq_ss00.txt')
axial_us = np.loadtxt(f'../../Gosha/Barions/Goshe/04axial_diq_us00.txt')
axial_ud = np.loadtxt(f'../../Gosha/Barions/Goshe/04axial_diq_ud00.txt')

scalar_ud = np.loadtxt(f'../../Gosha/Barions/Goshe/04scalar_diq_ud00.txt')
scalar_us = np.loadtxt(f'../../Gosha/Barions/Goshe/04scalar_diq_us00.txt')

m_u = np.loadtxt(f'../../Gosha/Barions/Goshe/04mueq_u_mass_phi_mu00.txt')
m_s = np.loadtxt(f'../../Gosha/Barions/Goshe/04mueq_smass_mu00.txt')

T_axial_ss = axial_ss[:300, 0]
M_axial_ss = axial_ss[:300, 1]

T_axial_us = axial_us[:300, 0]
M_axial_us = axial_us[:300, 1]

T_axial_ud = axial_ud[:300, 0]
M_axial_ud = axial_ud[:300, 1]

T_scalar_ud = scalar_ud[:300, 0]
M_scalar_ud = scalar_ud[:300, 1]

T_scalar_us = scalar_us[:300, 0]
M_scalar_us = scalar_us[:300, 1]

T_mu = m_u[:300, 0]
mq_u = m_u[:300, 2]

mq_s = m_s[:300, 2]

ax.plot(T_axial_ss, M_axial_ss, label='[ss]', linewidth=3)
ax.plot(T_axial_us, M_axial_us, label='[us]', linewidth=3, color='blue')
ax.plot(T_axial_ud, M_axial_ud, label='[ud]', linewidth=3, color='brown')
ax.plot(T_scalar_ud, M_scalar_ud, label='(ud)', linewidth=3, color='black')
ax.plot(T_scalar_us, M_scalar_us, label='(us)', linewidth=3, color='red')

ax.plot(T_mu, 2 * mq_u, '--', label='$2m_u', linewidth=2, color='black')
ax.plot(T_mu, mq_s + mq_u, '-.', label='$m_u + $m_s', linewidth=2, color='black')
ax.plot(T_mu, 2 * mq_s, ':', label='2$m_s', linewidth=2)
plt.text(0.03, 1.04, '[ss]', fontsize=18)
plt.text(0.03, 0.92, '[us]', fontsize=18)
plt.text(0.273, 0.83, '[ud]', fontsize=18)
plt.text(0.03, 0.60, '(ud)', fontsize=18)
plt.text(0.25, 0.64, '(us)', fontsize=18)

plt.text(0.23, 0.02, r'${2m_u}$', fontsize=18)
plt.text(0.235, 0.3, r'${m_u + m_s}$', fontsize=18)
plt.text(0.153, 1.1, r'${2m_s}$', fontsize=18)
ax.set(xlim=[0.003, 0.3], ylim=[0, 1.2])
ax.tick_params(which='major', direction='in', width=2, length=6)
for axis in ['top', 'bottom', 'left', 'right']:
    ax.spines[axis].set_linewidth(2)
plt.xlabel("Т, GeV", fontsize=18)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.ylabel("masses, GeV", fontsize=18)
plt.show()
# fig.savefig('Diquarks0.eps')
#
# proton = np.loadtxt(f'../../Gosha/Barions/Goshe/04proton00.txt')
# lambda_bar = np.loadtxt(f'../../Gosha/Barions/Goshe/04lambda_bar00.txt')
# psi = np.loadtxt(f'../../Gosha/Barions/Goshe/04psi_bar00.txt')
#
# T_p = proton[:, 0]
# M_p = proton[:182, 1]
# mu3 = proton[:, 2]
#
# T_l = lambda_bar[:, 0]
# l_b = lambda_bar[:197, 1]
# mq_u2_mq_s = lambda_bar[:, 2]
#
# T_psi = psi[:, 0]
# psi_bar = psi[:239, 1]
# mq_u2_mq_s_psi = psi[:, 2]
#
# ax.plot(T_p[:182], M_p, label='$p$', linewidth=3, color='blue')
# ax.plot(T_l[:197], l_b, label='$\Lambda^0$', linewidth=3, color='brown')
# ax.plot(T_psi[:239], psi_bar, label='$\Xi^-$', linewidth=3, color='green')
#
# ax.plot(T_p, mu3, '--', label='$2m_u$', linewidth=2, color='blue')
# ax.plot(T_psi, mq_u2_mq_s_psi, '-.', label='$2m_u$ + $m_s$', linewidth=2, color='green')
# ax.plot(T_l, mq_u2_mq_s, ':', label='$m_u$ + $2m_s$', linewidth=2, color='brown')
#
# plt.text(0.21, 1.2, '$m_u + 2m_s$', fontsize=18)
# plt.text(0.23, 0.3, '$2m_u + m_s$', fontsize=18)
# plt.text(0.21, 0.11, '$2m_u$ ', fontsize=18)
# plt.text(0.19, 0.9, '$\Lambda^0$', fontsize=20)
# plt.text(0.03, 1.29, '$\Xi^-$ ', fontsize=20,)
# plt.text(0.03, 0.92, '$p$ ', fontsize=20)
# ax.set(xlim=[0.003, 0.3], ylim=[0, 1.55])
# ax.tick_params( which='major', direction='in', width=2, length=6)
# for axis in ['top','bottom','left','right']:
#     ax.spines[axis].set_linewidth(2)
# plt.xlabel("Т, GeV", fontsize=18 )
# plt.xticks(fontsize=15)
# plt.yticks(fontsize=15)
# plt.ylabel("masses, GeV", fontsize=18)
# plt.show()
# fig.savefig('Barions0.eps')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

mu = ['00', '01', '03', '05', '07', '09', '11', '13', '15', '17', '19', '21', '23', '25', '26', '27', '28', '29',
      '30', '31', '32', '33', '34', '35']
plt.rcParams['text.latex.preamble'] = r"\usepackage{bm} \usepackage{amsmath} '\boldmath"
# '36', '37']
Tsc_list = []
fig, ax = plt.subplots()
# lambda_bar = np.loadtxt('../../../NJL/Goshe/Goshe/04axial_diq_ud00.txt')
# T_l = lambda_bar[:, 0]
# l_b = lambda_bar[:, 1]
# mq_u2_mq_s = lambda_bar[:, 2]
# plt.plot(T_l, l_b, label='line 1')
# plt.plot(T_l, mq_u2_mq_s, label='line 2')
# plt.show()


# for i in mu:
#     scalar = np.loadtxt(f'../../Gosha/Barions/Goshe/04axial_diq_ss{i}.txt')
#     # # if i in ['00', '01', '03', '05', '07',
#     #                                                                                   '09', '11'] else np.loadtxt(
#     #     f'../../Gosha/Barions/Goshe/04axial_diq{i}.txt')
#     m_u = np.loadtxt(f'../../Gosha/Barions/Goshe/04mueq_u_mass_phi_mu{i}.txt')
#     m_s = np.loadtxt(f'../../Gosha/Barions/Goshe/04mueq_smass_mu{i}.txt')
#     T_p = scalar[:, 0]
#     print(len(T_p))
#     m_ud_sc = scalar[0:280, 2] #if i in ['00', '01', '03', '05', '07', '09', '11'] else scalar[0:280, 2]
#     mq_u = 2 * m_u[0:280, 2]
#     mq_s = m_s[0:280, 2]
#     porog = mq_u + mq_s
#     Tsc_list.append(T_p[np.argmin(np.absolute(2*mq_s - porog))])
# # # print(Td_List)
# # print(T_lambda_bar_list)
#
# mu_2 = [0.00, 0.01, 0.03, 0.05, 0.07, 0.09, 0.11, 0.13, 0.15, 0.17, 0.19, 0.21, 0.23, 0.25, 0.26, 0.27, 0.28, 0.29,
#         0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37]
# # np.savetxt('T_mu_lambda.txt', np.c_[mu_2, T_lambda_bar_list], delimiter='    ')
# # np.savetxt('T_mu_protons.txt', np.c_[mu_2, Td_List], delimiter='    ')
# # np.savetxt('T_mu_psi.txt', np.c_[mu_2, T_psi_list], delimiter='    ')
#
# plt.scatter(mu, Tsc_list, c="g", label="axial_ss")
# # plt.scatter(mu_2, T_lambda_bar_list, c="r", label="lambda", marker="^")
# # plt.scatter(mu_2, T_psi_list, c="b", label="psi", marker="*")

#
axial_ss = np.loadtxt(f'../../Gosha/Barions/Goshe/04axial_diq_ss00.txt')
axial_us = np.loadtxt(f'../../Gosha/Barions/Goshe/04axial_diq_us00.txt')
axial_ud = np.loadtxt(f'../../Gosha/Barions/Goshe/04axial_diq_ud00.txt')

scalar_ud = np.loadtxt(f'../../Gosha/Barions/Goshe/04scalar_diq_ud00.txt')
scalar_us = np.loadtxt(f'../../Gosha/Barions/Goshe/04scalar_diq_us00.txt')

m_u = np.loadtxt(f'../../Gosha/Barions/Goshe/04mueq_u_mass_phi_mu00.txt')
m_s = np.loadtxt(f'../../Gosha/Barions/Goshe/04mueq_smass_mu00.txt')

T_axial_ss = axial_ss[:300, 0]
M_axial_ss = axial_ss[:300, 1]

T_axial_us = axial_us[:300, 0]
M_axial_us = axial_us[:300, 1]

T_axial_ud = axial_ud[:300, 0]
M_axial_ud = axial_ud[:300, 1]

T_scalar_ud = scalar_ud[:300, 0]
M_scalar_ud = scalar_ud[:300, 1]

T_scalar_us = scalar_us[:300, 0]
M_scalar_us = scalar_us[:300, 1]

T_mu = m_u[:300, 0]
mq_u = m_u[:300, 2]

mq_s = m_s[:300, 2]

ax.plot(T_axial_ss, M_axial_ss, label='[ss]', linewidth=3)
ax.plot(T_axial_us, M_axial_us, label='[us]', linewidth=3, color='blue')
ax.plot(T_axial_ud, M_axial_ud, label='[ud]', linewidth=3, color='brown')
ax.plot(T_scalar_ud, M_scalar_ud, label='(ud)', linewidth=3, color='black')
ax.plot(T_scalar_us, M_scalar_us, label='(us)', linewidth=3, color='red')

ax.plot(T_mu, 2 * mq_u, '--', label='$2m_u', linewidth=2, color='black')
ax.plot(T_mu, mq_s + mq_u, '-.', label='$m_u + $m_s', linewidth=2, color='black')
ax.plot(T_mu, 2 * mq_s, ':', label='2$m_s', linewidth=2)
plt.text(0.03, 1.04, '[ss]', fontsize=18)
plt.text(0.03, 0.92, '[us]', fontsize=18)
plt.text(0.273, 0.83, '[ud]', fontsize=18)
plt.text(0.03, 0.60, '(ud)', fontsize=18)
plt.text(0.25, 0.64, '(us)', fontsize=18)

plt.text(0.23, 0.02, r'${2m_u}$', fontsize=18)
plt.text(0.235, 0.3, r'${m_u + m_s}$', fontsize=18)
plt.text(0.153, 1.1, r'${2m_s}$', fontsize=18)
ax.set(xlim=[0.003, 0.3], ylim=[0, 1.2])
ax.tick_params(which='major', direction='in', width=2, length=6)
for axis in ['top', 'bottom', 'left', 'right']:
    ax.spines[axis].set_linewidth(2)
plt.xlabel("Т, GeV", fontsize=18)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.ylabel("masses, GeV", fontsize=18)
plt.show()
# fig.savefig('Diquarks0.eps')
#
# proton = np.loadtxt(f'../../Gosha/Barions/Goshe/04proton00.txt')
# lambda_bar = np.loadtxt(f'../../Gosha/Barions/Goshe/04lambda_bar00.txt')
# psi = np.loadtxt(f'../../Gosha/Barions/Goshe/04psi_bar00.txt')
#
# T_p = proton[:, 0]
# M_p = proton[:182, 1]
# mu3 = proton[:, 2]
#
# T_l = lambda_bar[:, 0]
# l_b = lambda_bar[:197, 1]
# mq_u2_mq_s = lambda_bar[:, 2]
#
# T_psi = psi[:, 0]
# psi_bar = psi[:239, 1]
# mq_u2_mq_s_psi = psi[:, 2]
#
# ax.plot(T_p[:182], M_p, label='$p$', linewidth=3, color='blue')
# ax.plot(T_l[:197], l_b, label='$\Lambda^0$', linewidth=3, color='brown')
# ax.plot(T_psi[:239], psi_bar, label='$\Xi^-$', linewidth=3, color='green')
#
# ax.plot(T_p, mu3, '--', label='$2m_u$', linewidth=2, color='blue')
# ax.plot(T_psi, mq_u2_mq_s_psi, '-.', label='$2m_u$ + $m_s$', linewidth=2, color='green')
# ax.plot(T_l, mq_u2_mq_s, ':', label='$m_u$ + $2m_s$', linewidth=2, color='brown')
#
# plt.text(0.21, 1.2, '$m_u + 2m_s$', fontsize=18)
# plt.text(0.23, 0.3, '$2m_u + m_s$', fontsize=18)
# plt.text(0.21, 0.11, '$2m_u$ ', fontsize=18)
# plt.text(0.19, 0.9, '$\Lambda^0$', fontsize=20)
# plt.text(0.03, 1.29, '$\Xi^-$ ', fontsize=20,)
# plt.text(0.03, 0.92, '$p$ ', fontsize=20)
# ax.set(xlim=[0.003, 0.3], ylim=[0, 1.55])
# ax.tick_params( which='major', direction='in', width=2, length=6)
# for axis in ['top','bottom','left','right']:
#     ax.spines[axis].set_linewidth(2)
# plt.xlabel("Т, GeV", fontsize=18 )
# plt.xticks(fontsize=15)
# plt.yticks(fontsize=15)
# plt.ylabel("masses, GeV", fontsize=18)
# plt.show()
# fig.savefig('Barions0.eps')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

mu = ['00', '01', '03', '05', '07', '09', '11', '13', '15', '17', '19', '21', '23', '25', '26', '27', '28', '29',
      '30', '31', '32', '33', '34', '35']
plt.rcParams['text.latex.preamble'] = r"\usepackage{bm} \usepackage{amsmath} '\boldmath"
# '36', '37']
Tsc_list = []
fig, ax = plt.subplots()
# lambda_bar = np.loadtxt('../../../NJL/Goshe/Goshe/04axial_diq_ud00.txt')
# T_l = lambda_bar[:, 0]
# l_b = lambda_bar[:, 1]
# mq_u2_mq_s = lambda_bar[:, 2]
# plt.plot(T_l, l_b, label='line 1')
# plt.plot(T_l, mq_u2_mq_s, label='line 2')
# plt.show()


# for i in mu:
#     scalar = np.loadtxt(f'../../Gosha/Barions/Goshe/04axial_diq_ss{i}.txt')
#     # # if i in ['00', '01', '03', '05', '07',
#     #                                                                                   '09', '11'] else np.loadtxt(
#     #     f'../../Gosha/Barions/Goshe/04axial_diq{i}.txt')
#     m_u = np.loadtxt(f'../../Gosha/Barions/Goshe/04mueq_u_mass_phi_mu{i}.txt')
#     m_s = np.loadtxt(f'../../Gosha/Barions/Goshe/04mueq_smass_mu{i}.txt')
#     T_p = scalar[:, 0]
#     print(len(T_p))
#     m_ud_sc = scalar[0:280, 2] #if i in ['00', '01', '03', '05', '07', '09', '11'] else scalar[0:280, 2]
#     mq_u = 2 * m_u[0:280, 2]
#     mq_s = m_s[0:280, 2]
#     porog = mq_u + mq_s
#     Tsc_list.append(T_p[np.argmin(np.absolute(2*mq_s - porog))])
# # # print(Td_List)
# # print(T_lambda_bar_list)
#
# mu_2 = [0.00, 0.01, 0.03, 0.05, 0.07, 0.09, 0.11, 0.13, 0.15, 0.17, 0.19, 0.21, 0.23, 0.25, 0.26, 0.27, 0.28, 0.29,
#         0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37]
# # np.savetxt('T_mu_lambda.txt', np.c_[mu_2, T_lambda_bar_list], delimiter='    ')
# # np.savetxt('T_mu_protons.txt', np.c_[mu_2, Td_List], delimiter='    ')
# # np.savetxt('T_mu_psi.txt', np.c_[mu_2, T_psi_list], delimiter='    ')
#
# plt.scatter(mu, Tsc_list, c="g", label="axial_ss")
# # plt.scatter(mu_2, T_lambda_bar_list, c="r", label="lambda", marker="^")
# # plt.scatter(mu_2, T_psi_list, c="b", label="psi", marker="*")

#
axial_ss = np.loadtxt(f'../../Gosha/Barions/Goshe/04axial_diq_ss00.txt')
axial_us = np.loadtxt(f'../../Gosha/Barions/Goshe/04axial_diq_us00.txt')
axial_ud = np.loadtxt(f'../../Gosha/Barions/Goshe/04axial_diq_ud00.txt')

scalar_ud = np.loadtxt(f'../../Gosha/Barions/Goshe/04scalar_diq_ud00.txt')
scalar_us = np.loadtxt(f'../../Gosha/Barions/Goshe/04scalar_diq_us00.txt')

m_u = np.loadtxt(f'../../Gosha/Barions/Goshe/04mueq_u_mass_phi_mu00.txt')
m_s = np.loadtxt(f'../../Gosha/Barions/Goshe/04mueq_smass_mu00.txt')

T_axial_ss = axial_ss[:300, 0]
M_axial_ss = axial_ss[:300, 1]

T_axial_us = axial_us[:300, 0]
M_axial_us = axial_us[:300, 1]

T_axial_ud = axial_ud[:300, 0]
M_axial_ud = axial_ud[:300, 1]

T_scalar_ud = scalar_ud[:300, 0]
M_scalar_ud = scalar_ud[:300, 1]

T_scalar_us = scalar_us[:300, 0]
M_scalar_us = scalar_us[:300, 1]

T_mu = m_u[:300, 0]
mq_u = m_u[:300, 2]

mq_s = m_s[:300, 2]

ax.plot(T_axial_ss, M_axial_ss, label='[ss]', linewidth=3)
ax.plot(T_axial_us, M_axial_us, label='[us]', linewidth=3, color='blue')
ax.plot(T_axial_ud, M_axial_ud, label='[ud]', linewidth=3, color='brown')
ax.plot(T_scalar_ud, M_scalar_ud, label='(ud)', linewidth=3, color='black')
ax.plot(T_scalar_us, M_scalar_us, label='(us)', linewidth=3, color='red')

ax.plot(T_mu, 2 * mq_u, '--', label='$2m_u', linewidth=2, color='black')
ax.plot(T_mu, mq_s + mq_u, '-.', label='$m_u + $m_s', linewidth=2, color='black')
ax.plot(T_mu, 2 * mq_s, ':', label='2$m_s', linewidth=2)
plt.text(0.03, 1.04, '[ss]', fontsize=18)
plt.text(0.03, 0.92, '[us]', fontsize=18)
plt.text(0.273, 0.83, '[ud]', fontsize=18)
plt.text(0.03, 0.60, '(ud)', fontsize=18)
plt.text(0.25, 0.64, '(us)', fontsize=18)

plt.text(0.23, 0.02, r'${2m_u}$', fontsize=18)
plt.text(0.235, 0.3, r'${m_u + m_s}$', fontsize=18)
plt.text(0.153, 1.1, r'${2m_s}$', fontsize=18)
ax.set(xlim=[0.003, 0.3], ylim=[0, 1.2])
ax.tick_params(which='major', direction='in', width=2, length=6)
for axis in ['top', 'bottom', 'left', 'right']:
    ax.spines[axis].set_linewidth(2)
plt.xlabel("Т, GeV", fontsize=18)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.ylabel("masses, GeV", fontsize=18)
plt.show()
# fig.savefig('Diquarks0.eps')
#
# proton = np.loadtxt(f'../../Gosha/Barions/Goshe/04proton00.txt')
# lambda_bar = np.loadtxt(f'../../Gosha/Barions/Goshe/04lambda_bar00.txt')
# psi = np.loadtxt(f'../../Gosha/Barions/Goshe/04psi_bar00.txt')
#
# T_p = proton[:, 0]
# M_p = proton[:182, 1]
# mu3 = proton[:, 2]
#
# T_l = lambda_bar[:, 0]
# l_b = lambda_bar[:197, 1]
# mq_u2_mq_s = lambda_bar[:, 2]
#
# T_psi = psi[:, 0]
# psi_bar = psi[:239, 1]
# mq_u2_mq_s_psi = psi[:, 2]
#
# ax.plot(T_p[:182], M_p, label='$p$', linewidth=3, color='blue')
# ax.plot(T_l[:197], l_b, label='$\Lambda^0$', linewidth=3, color='brown')
# ax.plot(T_psi[:239], psi_bar, label='$\Xi^-$', linewidth=3, color='green')
#
# ax.plot(T_p, mu3, '--', label='$2m_u$', linewidth=2, color='blue')
# ax.plot(T_psi, mq_u2_mq_s_psi, '-.', label='$2m_u$ + $m_s$', linewidth=2, color='green')
# ax.plot(T_l, mq_u2_mq_s, ':', label='$m_u$ + $2m_s$', linewidth=2, color='brown')
#
# plt.text(0.21, 1.2, '$m_u + 2m_s$', fontsize=18)
# plt.text(0.23, 0.3, '$2m_u + m_s$', fontsize=18)
# plt.text(0.21, 0.11, '$2m_u$ ', fontsize=18)
# plt.text(0.19, 0.9, '$\Lambda^0$', fontsize=20)
# plt.text(0.03, 1.29, '$\Xi^-$ ', fontsize=20,)
# plt.text(0.03, 0.92, '$p$ ', fontsize=20)
# ax.set(xlim=[0.003, 0.3], ylim=[0, 1.55])
# ax.tick_params( which='major', direction='in', width=2, length=6)
# for axis in ['top','bottom','left','right']:
#     ax.spines[axis].set_linewidth(2)
# plt.xlabel("Т, GeV", fontsize=18 )
# plt.xticks(fontsize=15)
# plt.yticks(fontsize=15)
# plt.ylabel("masses, GeV", fontsize=18)
# plt.show()
# fig.savefig('Barions0.eps')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

mu = ['00', '01', '03', '05', '07', '09', '11', '13', '15', '17', '19', '21', '23', '25', '26', '27', '28', '29',
      '30', '31', '32', '33', '34', '35']
plt.rcParams['text.latex.preamble'] = r"\usepackage{bm} \usepackage{amsmath} '\boldmath"
# '36', '37']
Tsc_list = []
fig, ax = plt.subplots()
# lambda_bar = np.loadtxt('../../../NJL/Goshe/Goshe/04axial_diq_ud00.txt')
# T_l = lambda_bar[:, 0]
# l_b = lambda_bar[:, 1]
# mq_u2_mq_s = lambda_bar[:, 2]
# plt.plot(T_l, l_b, label='line 1')
# plt.plot(T_l, mq_u2_mq_s, label='line 2')
# plt.show()


# for i in mu:
#     scalar = np.loadtxt(f'../../Gosha/Barions/Goshe/04axial_diq_ss{i}.txt')
#     # # if i in ['00', '01', '03', '05', '07',
#     #                                                                                   '09', '11'] else np.loadtxt(
#     #     f'../../Gosha/Barions/Goshe/04axial_diq{i}.txt')
#     m_u = np.loadtxt(f'../../Gosha/Barions/Goshe/04mueq_u_mass_phi_mu{i}.txt')
#     m_s = np.loadtxt(f'../../Gosha/Barions/Goshe/04mueq_smass_mu{i}.txt')
#     T_p = scalar[:, 0]
#     print(len(T_p))
#     m_ud_sc = scalar[0:280, 2] #if i in ['00', '01', '03', '05', '07', '09', '11'] else scalar[0:280, 2]
#     mq_u = 2 * m_u[0:280, 2]
#     mq_s = m_s[0:280, 2]
#     porog = mq_u + mq_s
#     Tsc_list.append(T_p[np.argmin(np.absolute(2*mq_s - porog))])
# # # print(Td_List)
# # print(T_lambda_bar_list)
#
# mu_2 = [0.00, 0.01, 0.03, 0.05, 0.07, 0.09, 0.11, 0.13, 0.15, 0.17, 0.19, 0.21, 0.23, 0.25, 0.26, 0.27, 0.28, 0.29,
#         0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37]
# # np.savetxt('T_mu_lambda.txt', np.c_[mu_2, T_lambda_bar_list], delimiter='    ')
# # np.savetxt('T_mu_protons.txt', np.c_[mu_2, Td_List], delimiter='    ')
# # np.savetxt('T_mu_psi.txt', np.c_[mu_2, T_psi_list], delimiter='    ')
#
# plt.scatter(mu, Tsc_list, c="g", label="axial_ss")
# # plt.scatter(mu_2, T_lambda_bar_list, c="r", label="lambda", marker="^")
# # plt.scatter(mu_2, T_psi_list, c="b", label="psi", marker="*")

#
axial_ss = np.loadtxt(f'../../Gosha/Barions/Goshe/04axial_diq_ss00.txt')
axial_us = np.loadtxt(f'../../Gosha/Barions/Goshe/04axial_diq_us00.txt')
axial_ud = np.loadtxt(f'../../Gosha/Barions/Goshe/04axial_diq_ud00.txt')

scalar_ud = np.loadtxt(f'../../Gosha/Barions/Goshe/04scalar_diq_ud00.txt')
scalar_us = np.loadtxt(f'../../Gosha/Barions/Goshe/04scalar_diq_us00.txt')

m_u = np.loadtxt(f'../../Gosha/Barions/Goshe/04mueq_u_mass_phi_mu00.txt')
m_s = np.loadtxt(f'../../Gosha/Barions/Goshe/04mueq_smass_mu00.txt')

T_axial_ss = axial_ss[:300, 0]
M_axial_ss = axial_ss[:300, 1]

T_axial_us = axial_us[:300, 0]
M_axial_us = axial_us[:300, 1]

T_axial_ud = axial_ud[:300, 0]
M_axial_ud = axial_ud[:300, 1]

T_scalar_ud = scalar_ud[:300, 0]
M_scalar_ud = scalar_ud[:300, 1]

T_scalar_us = scalar_us[:300, 0]
M_scalar_us = scalar_us[:300, 1]

T_mu = m_u[:300, 0]
mq_u = m_u[:300, 2]

mq_s = m_s[:300, 2]

ax.plot(T_axial_ss, M_axial_ss, label='[ss]', linewidth=3)
ax.plot(T_axial_us, M_axial_us, label='[us]', linewidth=3, color='blue')
ax.plot(T_axial_ud, M_axial_ud, label='[ud]', linewidth=3, color='brown')
ax.plot(T_scalar_ud, M_scalar_ud, label='(ud)', linewidth=3, color='black')
ax.plot(T_scalar_us, M_scalar_us, label='(us)', linewidth=3, color='red')

ax.plot(T_mu, 2 * mq_u, '--', label='$2m_u', linewidth=2, color='black')
ax.plot(T_mu, mq_s + mq_u, '-.', label='$m_u + $m_s', linewidth=2, color='black')
ax.plot(T_mu, 2 * mq_s, ':', label='2$m_s', linewidth=2)
plt.text(0.03, 1.04, '[ss]', fontsize=18)
plt.text(0.03, 0.92, '[us]', fontsize=18)
plt.text(0.273, 0.83, '[ud]', fontsize=18)
plt.text(0.03, 0.60, '(ud)', fontsize=18)
plt.text(0.25, 0.64, '(us)', fontsize=18)

plt.text(0.23, 0.02, r'${2m_u}$', fontsize=18)
plt.text(0.235, 0.3, r'${m_u + m_s}$', fontsize=18)
plt.text(0.153, 1.1, r'${2m_s}$', fontsize=18)
ax.set(xlim=[0.003, 0.3], ylim=[0, 1.2])
ax.tick_params(which='major', direction='in', width=2, length=6)
for axis in ['top', 'bottom', 'left', 'right']:
    ax.spines[axis].set_linewidth(2)
plt.xlabel("Т, GeV", fontsize=18)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.ylabel("masses, GeV", fontsize=18)
plt.show()
# fig.savefig('Diquarks0.eps')
#
# proton = np.loadtxt(f'../../Gosha/Barions/Goshe/04proton00.txt')
# lambda_bar = np.loadtxt(f'../../Gosha/Barions/Goshe/04lambda_bar00.txt')
# psi = np.loadtxt(f'../../Gosha/Barions/Goshe/04psi_bar00.txt')
#
# T_p = proton[:, 0]
# M_p = proton[:182, 1]
# mu3 = proton[:, 2]
#
# T_l = lambda_bar[:, 0]
# l_b = lambda_bar[:197, 1]
# mq_u2_mq_s = lambda_bar[:, 2]
#
# T_psi = psi[:, 0]
# psi_bar = psi[:239, 1]
# mq_u2_mq_s_psi = psi[:, 2]
#
# ax.plot(T_p[:182], M_p, label='$p$', linewidth=3, color='blue')
# ax.plot(T_l[:197], l_b, label='$\Lambda^0$', linewidth=3, color='brown')
# ax.plot(T_psi[:239], psi_bar, label='$\Xi^-$', linewidth=3, color='green')
#
# ax.plot(T_p, mu3, '--', label='$2m_u$', linewidth=2, color='blue')
# ax.plot(T_psi, mq_u2_mq_s_psi, '-.', label='$2m_u$ + $m_s$', linewidth=2, color='green')
# ax.plot(T_l, mq_u2_mq_s, ':', label='$m_u$ + $2m_s$', linewidth=2, color='brown')
#
# plt.text(0.21, 1.2, '$m_u + 2m_s$', fontsize=18)
# plt.text(0.23, 0.3, '$2m_u + m_s$', fontsize=18)
# plt.text(0.21, 0.11, '$2m_u$ ', fontsize=18)
# plt.text(0.19, 0.9, '$\Lambda^0$', fontsize=20)
# plt.text(0.03, 1.29, '$\Xi^-$ ', fontsize=20,)
# plt.text(0.03, 0.92, '$p$ ', fontsize=20)
# ax.set(xlim=[0.003, 0.3], ylim=[0, 1.55])
# ax.tick_params( which='major', direction='in', width=2, length=6)
# for axis in ['top','bottom','left','right']:
#     ax.spines[axis].set_linewidth(2)
# plt.xlabel("Т, GeV", fontsize=18 )
# plt.xticks(fontsize=15)
# plt.yticks(fontsize=15)
# plt.ylabel("masses, GeV", fontsize=18)
# plt.show()
# fig.savefig('Barions0.eps')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

mu = ['00', '01', '03', '05', '07', '09', '11', '13', '15', '17', '19', '21', '23', '25', '26', '27', '28', '29',
      '30', '31', '32', '33', '34', '35']
plt.rcParams['text.latex.preamble'] = r"\usepackage{bm} \usepackage{amsmath} '\boldmath"
# '36', '37']
Tsc_list = []
fig, ax = plt.subplots()
# lambda_bar = np.loadtxt('../../../NJL/Goshe/Goshe/04axial_diq_ud00.txt')
# T_l = lambda_bar[:, 0]
# l_b = lambda_bar[:, 1]
# mq_u2_mq_s = lambda_bar[:, 2]
# plt.plot(T_l, l_b, label='line 1')
# plt.plot(T_l, mq_u2_mq_s, label='line 2')
# plt.show()


# for i in mu:
#     scalar = np.loadtxt(f'../../Gosha/Barions/Goshe/04axial_diq_ss{i}.txt')
#     # # if i in ['00', '01', '03', '05', '07',
#     #                                                                                   '09', '11'] else np.loadtxt(
#     #     f'../../Gosha/Barions/Goshe/04axial_diq{i}.txt')
#     m_u = np.loadtxt(f'../../Gosha/Barions/Goshe/04mueq_u_mass_phi_mu{i}.txt')
#     m_s = np.loadtxt(f'../../Gosha/Barions/Goshe/04mueq_smass_mu{i}.txt')
#     T_p = scalar[:, 0]
#     print(len(T_p))
#     m_ud_sc = scalar[0:280, 2] #if i in ['00', '01', '03', '05', '07', '09', '11'] else scalar[0:280, 2]
#     mq_u = 2 * m_u[0:280, 2]
#     mq_s = m_s[0:280, 2]
#     porog = mq_u + mq_s
#     Tsc_list.append(T_p[np.argmin(np.absolute(2*mq_s - porog))])
# # # print(Td_List)
# # print(T_lambda_bar_list)
#
# mu_2 = [0.00, 0.01, 0.03, 0.05, 0.07, 0.09, 0.11, 0.13, 0.15, 0.17, 0.19, 0.21, 0.23, 0.25, 0.26, 0.27, 0.28, 0.29,
#         0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37]
# # np.savetxt('T_mu_lambda.txt', np.c_[mu_2, T_lambda_bar_list], delimiter='    ')
# # np.savetxt('T_mu_protons.txt', np.c_[mu_2, Td_List], delimiter='    ')
# # np.savetxt('T_mu_psi.txt', np.c_[mu_2, T_psi_list], delimiter='    ')
#
# plt.scatter(mu, Tsc_list, c="g", label="axial_ss")
# # plt.scatter(mu_2, T_lambda_bar_list, c="r", label="lambda", marker="^")
# # plt.scatter(mu_2, T_psi_list, c="b", label="psi", marker="*")

#
axial_ss = np.loadtxt(f'../../Gosha/Barions/Goshe/04axial_diq_ss00.txt')
axial_us = np.loadtxt(f'../../Gosha/Barions/Goshe/04axial_diq_us00.txt')
axial_ud = np.loadtxt(f'../../Gosha/Barions/Goshe/04axial_diq_ud00.txt')

scalar_ud = np.loadtxt(f'../../Gosha/Barions/Goshe/04scalar_diq_ud00.txt')
scalar_us = np.loadtxt(f'../../Gosha/Barions/Goshe/04scalar_diq_us00.txt')

m_u = np.loadtxt(f'../../Gosha/Barions/Goshe/04mueq_u_mass_phi_mu00.txt')
m_s = np.loadtxt(f'../../Gosha/Barions/Goshe/04mueq_smass_mu00.txt')

T_axial_ss = axial_ss[:300, 0]
M_axial_ss = axial_ss[:300, 1]

T_axial_us = axial_us[:300, 0]
M_axial_us = axial_us[:300, 1]

T_axial_ud = axial_ud[:300, 0]
M_axial_ud = axial_ud[:300, 1]

T_scalar_ud = scalar_ud[:300, 0]
M_scalar_ud = scalar_ud[:300, 1]

T_scalar_us = scalar_us[:300, 0]
M_scalar_us = scalar_us[:300, 1]

T_mu = m_u[:300, 0]
mq_u = m_u[:300, 2]

mq_s = m_s[:300, 2]

ax.plot(T_axial_ss, M_axial_ss, label='[ss]', linewidth=3)
ax.plot(T_axial_us, M_axial_us, label='[us]', linewidth=3, color='blue')
ax.plot(T_axial_ud, M_axial_ud, label='[ud]', linewidth=3, color='brown')
ax.plot(T_scalar_ud, M_scalar_ud, label='(ud)', linewidth=3, color='black')
ax.plot(T_scalar_us, M_scalar_us, label='(us)', linewidth=3, color='red')

ax.plot(T_mu, 2 * mq_u, '--', label='$2m_u', linewidth=2, color='black')
ax.plot(T_mu, mq_s + mq_u, '-.', label='$m_u + $m_s', linewidth=2, color='black')
ax.plot(T_mu, 2 * mq_s, ':', label='2$m_s', linewidth=2)
plt.text(0.03, 1.04, '[ss]', fontsize=18)
plt.text(0.03, 0.92, '[us]', fontsize=18)
plt.text(0.273, 0.83, '[ud]', fontsize=18)
plt.text(0.03, 0.60, '(ud)', fontsize=18)
plt.text(0.25, 0.64, '(us)', fontsize=18)

plt.text(0.23, 0.02, r'${2m_u}$', fontsize=18)
plt.text(0.235, 0.3, r'${m_u + m_s}$', fontsize=18)
plt.text(0.153, 1.1, r'${2m_s}$', fontsize=18)
ax.set(xlim=[0.003, 0.3], ylim=[0, 1.2])
ax.tick_params(which='major', direction='in', width=2, length=6)
for axis in ['top', 'bottom', 'left', 'right']:
    ax.spines[axis].set_linewidth(2)
plt.xlabel("Т, GeV", fontsize=18)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.ylabel("masses, GeV", fontsize=18)
plt.show()
# fig.savefig('Diquarks0.eps')
#
# proton = np.loadtxt(f'../../Gosha/Barions/Goshe/04proton00.txt')
# lambda_bar = np.loadtxt(f'../../Gosha/Barions/Goshe/04lambda_bar00.txt')
# psi = np.loadtxt(f'../../Gosha/Barions/Goshe/04psi_bar00.txt')
#
# T_p = proton[:, 0]
# M_p = proton[:182, 1]
# mu3 = proton[:, 2]
#
# T_l = lambda_bar[:, 0]
# l_b = lambda_bar[:197, 1]
# mq_u2_mq_s = lambda_bar[:, 2]
#
# T_psi = psi[:, 0]
# psi_bar = psi[:239, 1]
# mq_u2_mq_s_psi = psi[:, 2]
#
# ax.plot(T_p[:182], M_p, label='$p$', linewidth=3, color='blue')
# ax.plot(T_l[:197], l_b, label='$\Lambda^0$', linewidth=3, color='brown')
# ax.plot(T_psi[:239], psi_bar, label='$\Xi^-$', linewidth=3, color='green')
#
# ax.plot(T_p, mu3, '--', label='$2m_u$', linewidth=2, color='blue')
# ax.plot(T_psi, mq_u2_mq_s_psi, '-.', label='$2m_u$ + $m_s$', linewidth=2, color='green')
# ax.plot(T_l, mq_u2_mq_s, ':', label='$m_u$ + $2m_s$', linewidth=2, color='brown')
#
# plt.text(0.21, 1.2, '$m_u + 2m_s$', fontsize=18)
# plt.text(0.23, 0.3, '$2m_u + m_s$', fontsize=18)
# plt.text(0.21, 0.11, '$2m_u$ ', fontsize=18)
# plt.text(0.19, 0.9, '$\Lambda^0$', fontsize=20)
# plt.text(0.03, 1.29, '$\Xi^-$ ', fontsize=20,)
# plt.text(0.03, 0.92, '$p$ ', fontsize=20)
# ax.set(xlim=[0.003, 0.3], ylim=[0, 1.55])
# ax.tick_params( which='major', direction='in', width=2, length=6)
# for axis in ['top','bottom','left','right']:
#     ax.spines[axis].set_linewidth(2)
# plt.xlabel("Т, GeV", fontsize=18 )
# plt.xticks(fontsize=15)
# plt.yticks(fontsize=15)
# plt.ylabel("masses, GeV", fontsize=18)
# plt.show()
# fig.savefig('Barions0.eps')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

mu = ['00', '01', '03', '05', '07', '09', '11', '13', '15', '17', '19', '21', '23', '25', '26', '27', '28', '29',
      '30', '31', '32', '33', '34', '35']
plt.rcParams['text.latex.preamble'] = r"\usepackage{bm} \usepackage{amsmath} '\boldmath"
# '36', '37']
Tsc_list = []
fig, ax = plt.subplots()
# lambda_bar = np.loadtxt('../../../NJL/Goshe/Goshe/04axial_diq_ud00.txt')
# T_l = lambda_bar[:, 0]
# l_b = lambda_bar[:, 1]
# mq_u2_mq_s = lambda_bar[:, 2]
# plt.plot(T_l, l_b, label='line 1')
# plt.plot(T_l, mq_u2_mq_s, label='line 2')
# plt.show()


# for i in mu:
#     scalar = np.loadtxt(f'../../Gosha/Barions/Goshe/04axial_diq_ss{i}.txt')
#     # # if i in ['00', '01', '03', '05', '07',
#     #                                                                                   '09', '11'] else np.loadtxt(
#     #     f'../../Gosha/Barions/Goshe/04axial_diq{i}.txt')
#     m_u = np.loadtxt(f'../../Gosha/Barions/Goshe/04mueq_u_mass_phi_mu{i}.txt')
#     m_s = np.loadtxt(f'../../Gosha/Barions/Goshe/04mueq_smass_mu{i}.txt')
#     T_p = scalar[:, 0]
#     print(len(T_p))
#     m_ud_sc = scalar[0:280, 2] #if i in ['00', '01', '03', '05', '07', '09', '11'] else scalar[0:280, 2]
#     mq_u = 2 * m_u[0:280, 2]
#     mq_s = m_s[0:280, 2]
#     porog = mq_u + mq_s
#     Tsc_list.append(T_p[np.argmin(np.absolute(2*mq_s - porog))])
# # # print(Td_List)
# # print(T_lambda_bar_list)
#
# mu_2 = [0.00, 0.01, 0.03, 0.05, 0.07, 0.09, 0.11, 0.13, 0.15, 0.17, 0.19, 0.21, 0.23, 0.25, 0.26, 0.27, 0.28, 0.29,
#         0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37]
# # np.savetxt('T_mu_lambda.txt', np.c_[mu_2, T_lambda_bar_list], delimiter='    ')
# # np.savetxt('T_mu_protons.txt', np.c_[mu_2, Td_List], delimiter='    ')
# # np.savetxt('T_mu_psi.txt', np.c_[mu_2, T_psi_list], delimiter='    ')
#
# plt.scatter(mu, Tsc_list, c="g", label="axial_ss")
# # plt.scatter(mu_2, T_lambda_bar_list, c="r", label="lambda", marker="^")
# # plt.scatter(mu_2, T_psi_list, c="b", label="psi", marker="*")

#
axial_ss = np.loadtxt(f'../../Gosha/Barions/Goshe/04axial_diq_ss00.txt')
axial_us = np.loadtxt(f'../../Gosha/Barions/Goshe/04axial_diq_us00.txt')
axial_ud = np.loadtxt(f'../../Gosha/Barions/Goshe/04axial_diq_ud00.txt')

scalar_ud = np.loadtxt(f'../../Gosha/Barions/Goshe/04scalar_diq_ud00.txt')
scalar_us = np.loadtxt(f'../../Gosha/Barions/Goshe/04scalar_diq_us00.txt')

m_u = np.loadtxt(f'../../Gosha/Barions/Goshe/04mueq_u_mass_phi_mu00.txt')
m_s = np.loadtxt(f'../../Gosha/Barions/Goshe/04mueq_smass_mu00.txt')

T_axial_ss = axial_ss[:300, 0]
M_axial_ss = axial_ss[:300, 1]

T_axial_us = axial_us[:300, 0]
M_axial_us = axial_us[:300, 1]

T_axial_ud = axial_ud[:300, 0]
M_axial_ud = axial_ud[:300, 1]

T_scalar_ud = scalar_ud[:300, 0]
M_scalar_ud = scalar_ud[:300, 1]

T_scalar_us = scalar_us[:300, 0]
M_scalar_us = scalar_us[:300, 1]

T_mu = m_u[:300, 0]
mq_u = m_u[:300, 2]

mq_s = m_s[:300, 2]

ax.plot(T_axial_ss, M_axial_ss, label='[ss]', linewidth=3)
ax.plot(T_axial_us, M_axial_us, label='[us]', linewidth=3, color='blue')
ax.plot(T_axial_ud, M_axial_ud, label='[ud]', linewidth=3, color='brown')
ax.plot(T_scalar_ud, M_scalar_ud, label='(ud)', linewidth=3, color='black')
ax.plot(T_scalar_us, M_scalar_us, label='(us)', linewidth=3, color='red')

ax.plot(T_mu, 2 * mq_u, '--', label='$2m_u', linewidth=2, color='black')
ax.plot(T_mu, mq_s + mq_u, '-.', label='$m_u + $m_s', linewidth=2, color='black')
ax.plot(T_mu, 2 * mq_s, ':', label='2$m_s', linewidth=2)
plt.text(0.03, 1.04, '[ss]', fontsize=18)
plt.text(0.03, 0.92, '[us]', fontsize=18)
plt.text(0.273, 0.83, '[ud]', fontsize=18)
plt.text(0.03, 0.60, '(ud)', fontsize=18)
plt.text(0.25, 0.64, '(us)', fontsize=18)

plt.text(0.23, 0.02, r'${2m_u}$', fontsize=18)
plt.text(0.235, 0.3, r'${m_u + m_s}$', fontsize=18)
plt.text(0.153, 1.1, r'${2m_s}$', fontsize=18)
ax.set(xlim=[0.003, 0.3], ylim=[0, 1.2])
ax.tick_params(which='major', direction='in', width=2, length=6)
for axis in ['top', 'bottom', 'left', 'right']:
    ax.spines[axis].set_linewidth(2)
plt.xlabel("Т, GeV", fontsize=18)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.ylabel("masses, GeV", fontsize=18)
plt.show()
# fig.savefig('Diquarks0.eps')
#
# proton = np.loadtxt(f'../../Gosha/Barions/Goshe/04proton00.txt')
# lambda_bar = np.loadtxt(f'../../Gosha/Barions/Goshe/04lambda_bar00.txt')
# psi = np.loadtxt(f'../../Gosha/Barions/Goshe/04psi_bar00.txt')
#
# T_p = proton[:, 0]
# M_p = proton[:182, 1]
# mu3 = proton[:, 2]
#
# T_l = lambda_bar[:, 0]
# l_b = lambda_bar[:197, 1]
# mq_u2_mq_s = lambda_bar[:, 2]
#
# T_psi = psi[:, 0]
# psi_bar = psi[:239, 1]
# mq_u2_mq_s_psi = psi[:, 2]
#
# ax.plot(T_p[:182], M_p, label='$p$', linewidth=3, color='blue')
# ax.plot(T_l[:197], l_b, label='$\Lambda^0$', linewidth=3, color='brown')
# ax.plot(T_psi[:239], psi_bar, label='$\Xi^-$', linewidth=3, color='green')
#
# ax.plot(T_p, mu3, '--', label='$2m_u$', linewidth=2, color='blue')
# ax.plot(T_psi, mq_u2_mq_s_psi, '-.', label='$2m_u$ + $m_s$', linewidth=2, color='green')
# ax.plot(T_l, mq_u2_mq_s, ':', label='$m_u$ + $2m_s$', linewidth=2, color='brown')
#
# plt.text(0.21, 1.2, '$m_u + 2m_s$', fontsize=18)
# plt.text(0.23, 0.3, '$2m_u + m_s$', fontsize=18)
# plt.text(0.21, 0.11, '$2m_u$ ', fontsize=18)
# plt.text(0.19, 0.9, '$\Lambda^0$', fontsize=20)
# plt.text(0.03, 1.29, '$\Xi^-$ ', fontsize=20,)
# plt.text(0.03, 0.92, '$p$ ', fontsize=20)
# ax.set(xlim=[0.003, 0.3], ylim=[0, 1.55])
# ax.tick_params( which='major', direction='in', width=2, length=6)
# for axis in ['top','bottom','left','right']:
#     ax.spines[axis].set_linewidth(2)
# plt.xlabel("Т, GeV", fontsize=18 )
# plt.xticks(fontsize=15)
# plt.yticks(fontsize=15)
# plt.ylabel("masses, GeV", fontsize=18)
# plt.show()
# fig.savefig('Barions0.eps')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

mu = ['00', '01', '03', '05', '07', '09', '11', '13', '15', '17', '19', '21', '23', '25', '26', '27', '28', '29',
      '30', '31', '32', '33', '34', '35']
plt.rcParams['text.latex.preamble'] = r"\usepackage{bm} \usepackage{amsmath} '\boldmath"
# '36', '37']
Tsc_list = []
fig, ax = plt.subplots()
# lambda_bar = np.loadtxt('../../../NJL/Goshe/Goshe/04axial_diq_ud00.txt')
# T_l = lambda_bar[:, 0]
# l_b = lambda_bar[:, 1]
# mq_u2_mq_s = lambda_bar[:, 2]
# plt.plot(T_l, l_b, label='line 1')
# plt.plot(T_l, mq_u2_mq_s, label='line 2')
# plt.show()


# for i in mu:
#     scalar = np.loadtxt(f'../../Gosha/Barions/Goshe/04axial_diq_ss{i}.txt')
#     # # if i in ['00', '01', '03', '05', '07',
#     #                                                                                   '09', '11'] else np.loadtxt(
#     #     f'../../Gosha/Barions/Goshe/04axial_diq{i}.txt')
#     m_u = np.loadtxt(f'../../Gosha/Barions/Goshe/04mueq_u_mass_phi_mu{i}.txt')
#     m_s = np.loadtxt(f'../../Gosha/Barions/Goshe/04mueq_smass_mu{i}.txt')
#     T_p = scalar[:, 0]
#     print(len(T_p))
#     m_ud_sc = scalar[0:280, 2] #if i in ['00', '01', '03', '05', '07', '09', '11'] else scalar[0:280, 2]
#     mq_u = 2 * m_u[0:280, 2]
#     mq_s = m_s[0:280, 2]
#     porog = mq_u + mq_s
#     Tsc_list.append(T_p[np.argmin(np.absolute(2*mq_s - porog))])
# # # print(Td_List)
# # print(T_lambda_bar_list)
#
# mu_2 = [0.00, 0.01, 0.03, 0.05, 0.07, 0.09, 0.11, 0.13, 0.15, 0.17, 0.19, 0.21, 0.23, 0.25, 0.26, 0.27, 0.28, 0.29,
#         0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37]
# # np.savetxt('T_mu_lambda.txt', np.c_[mu_2, T_lambda_bar_list], delimiter='    ')
# # np.savetxt('T_mu_protons.txt', np.c_[mu_2, Td_List], delimiter='    ')
# # np.savetxt('T_mu_psi.txt', np.c_[mu_2, T_psi_list], delimiter='    ')
#
# plt.scatter(mu, Tsc_list, c="g", label="axial_ss")
# # plt.scatter(mu_2, T_lambda_bar_list, c="r", label="lambda", marker="^")
# # plt.scatter(mu_2, T_psi_list, c="b", label="psi", marker="*")

#
# axial_ss = np.loadtxt(f'../../Gosha/Barions/Goshe/04axial_diq_ss00.txt')
# axial_us = np.loadtxt(f'../../Gosha/Barions/Goshe/04axial_diq_us00.txt')
# axial_ud = np.loadtxt(f'../../Gosha/Barions/Goshe/04axial_diq_ud00.txt')
#
# scalar_ud = np.loadtxt(f'../../Gosha/Barions/Goshe/04scalar_diq_ud00.txt')
# scalar_us = np.loadtxt(f'../../Gosha/Barions/Goshe/04scalar_diq_us00.txt')
#
# m_u = np.loadtxt(f'../../Gosha/Barions/Goshe/04mueq_u_mass_phi_mu00.txt')
# m_s = np.loadtxt(f'../../Gosha/Barions/Goshe/04mueq_smass_mu00.txt')
#
# T_axial_ss = axial_ss[:300, 0]
# M_axial_ss = axial_ss[:300, 1]
#
# T_axial_us = axial_us[:300, 0]
# M_axial_us = axial_us[:300, 1]
#
# T_axial_ud = axial_ud[:300, 0]
# M_axial_ud = axial_ud[:300, 1]
#
# T_scalar_ud = scalar_ud[:300, 0]
# M_scalar_ud = scalar_ud[:300, 1]
#
# T_scalar_us = scalar_us[:300, 0]
# M_scalar_us = scalar_us[:300, 1]
#
# T_mu = m_u[:300, 0]
# mq_u = m_u[:300, 2]
#
# mq_s = m_s[:300, 2]
#
# ax.plot(T_axial_ss, M_axial_ss, label='[ss]', linewidth=3, color='green')
# ax.plot(T_axial_us, M_axial_us, label='[us]', linewidth=3, color='blue')
# ax.plot(T_axial_ud, M_axial_ud, label='[ud]', linewidth=3, color='brown')
# ax.plot(T_scalar_ud, M_scalar_ud, label='(ud)', linewidth=3, color='black')
# ax.plot(T_scalar_us, M_scalar_us, label='(us)', linewidth=3, color='red')
#
# ax.plot(T_mu, 2 * mq_u, '--', label='$2m_u', linewidth=2, color='black')
# ax.plot(T_mu, mq_s + mq_u, '-.', label='$m_u + $m_s', linewidth=2, color='red')
# ax.plot(T_mu, 2 * mq_s, ':', label='2$m_s', linewidth=2, color='green')
# plt.text(0.03, 1.1, '[ss]', fontsize=18)
# plt.text(0.03, 0.92, '[us]', fontsize=18)
# plt.text(0.23, 0.67, '[ud]', fontsize=18)
# plt.text(0.03, 0.60, '(ud)', fontsize=18)
# plt.text(0.03, 0.75, '(us)', fontsize=18)
#
# plt.text(0.23, 0.02, r'${2m_u}$', fontsize=18)
# plt.text(0.235, 0.3, r'${m_u + m_s}$', fontsize=18)
# plt.text(0.153, 1.1, r'${2m_s}$', fontsize=18)
# ax.set(xlim=[0.003, 0.3], ylim=[0, 1.2])
# ax.tick_params(which='major', direction='in', width=2, length=6)
# for axis in ['top', 'bottom', 'left', 'right']:
#     ax.spines[axis].set_linewidth(2)
# plt.xlabel("Т, GeV", fontsize=18)
# plt.xticks(fontsize=15)
# plt.yticks(fontsize=15)
# plt.ylabel("masses, GeV", fontsize=18)
# plt.show()
# fig.savefig('Diquarks0.eps')
#
proton = np.loadtxt(f'../../Gosha/Barions/Goshe/04proton00.txt')
lambda_bar = np.loadtxt(f'../../Gosha/Barions/Goshe/04lambda_bar00.txt')
psi = np.loadtxt(f'../../Gosha/Barions/Goshe/04psi_bar00.txt')

T_p = proton[:, 0]
M_p = proton[:182, 1]
mu3 = proton[:, 2]

T_l = lambda_bar[:, 0]
l_b = lambda_bar[:197, 1]
mq_u2_mq_s = lambda_bar[:, 2]

T_psi = psi[:, 0]
psi_bar = psi[:239, 1]
mq_u2_mq_s_psi = psi[:, 2]

ax.plot(T_p[:182], M_p, label='$p$', linewidth=3, color='blue')
ax.plot(T_l[:197], l_b, label='$\Lambda^0$', linewidth=3, color='brown')
ax.plot(T_psi[:239], psi_bar, label='$\Xi^-$', linewidth=3, color='green')

ax.plot(T_p, mu3,  , label='$2m_u$', linewidth=2, color='blue')
ax.plot(T_psi, mq_u2_mq_s_psi, '-.', label='$2m_u$ + $m_s$', linewidth=2, color='green')
ax.plot(T_l, mq_u2_mq_s, ':', label='$m_u$ + $2m_s$', linewidth=2, color='brown')

plt.text(0.21, 1.2, '$m_u + 2m_s$', fontsize=18)
plt.text(0.23, 0.3, '$2m_u + m_s$', fontsize=18)
plt.text(0.21, 0.11, '$2m_u$ ', fontsize=18)
plt.text(0.19, 0.9, '$\Lambda^0$', fontsize=20)
plt.text(0.03, 1.29, '$\Xi^-$ ', fontsize=20,)
plt.text(0.03, 0.92, '$p$ ', fontsize=20)
ax.set(xlim=[0.003, 0.3], ylim=[0, 1.55])
ax.tick_params( which='major', direction='in', width=2, length=6)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(2)
plt.xlabel("Т, GeV", fontsize=18 )
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.ylabel("masses, GeV", fontsize=18)
plt.show()
fig.savefig('Barions0.eps')
