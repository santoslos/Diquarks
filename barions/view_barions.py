import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

fig, ax = plt.subplots()
mu = ['00', '01', '03', '05', '07', '09', '11',
      '13', '15', '17', '19', '21', '23', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35']

# q = np.loadtxt(f'../../Gosha/Barions/Goshe/Mu/04mueq_dmass_mu01.txt')
# ax.plot(q[:, 3], q[:, 2])
# plt.show()
#
# Td_List = []
# T_lambda_bar_list = []
# T_psi_list = []
# mueq_u_eb = np.loadtxt(f'../../Gosha/Barions/Goshe/Mu/newParametrs/04mueq_u_mass_phi_mu00_EB.txt')
# mueq_d_eb = np.loadtxt(f'../../Gosha/Barions/Goshe/Mu/newParametrs/04mueq_dmass_mu00_EB.txt')
# mueq_s_eb = np.loadtxt(f'../../Gosha/Barions/Goshe/Mu/newParametrs/04mueq_smass_mu00_EB.txt')
# mueq_u_p1 = np.loadtxt(f'../../Gosha/Barions/Goshe/Mu/newParametrs/04mueq_u_mass_phi_mu00.txt')
# mueq_d_p1 = np.loadtxt(f'../../Gosha/Barions/Goshe/Mu/newParametrs/04mueq_dmass_mu00.txt')
# mueq_s_1 = np.loadtxt(f'../../Gosha/Barions/Goshe/Mu/newParametrs/04mueq_smass_mu00.txt')
# u_q= np.loadtxt(f'../../Gosha/Barions/Goshe/04mueq_u_mass_phi_mu37.txt')
# s_q= np.loadtxt(f'../../Gosha/Barions/Goshe/04mueq_smass_mu37.txt')
# d_q= np.loadtxt(f'../../Gosha/Barions/Goshe/04mueq_dmass_mu37.txt')
# mdiq_ud_A = np.loadtxt(f'../../Gosha/Barions/Goshe/Mu/04scalar_diq_ud_T01_A.txt')
# mdiq_ud_B = np.loadtxt(f'../../Gosha/Barions/Goshe/Mu/04scalar_diq_ud_T01_B.txt')
# mdiq_ud = np.loadtxt(f'../../Gosha/Barions/Goshe/Mu/04scalar_diq_ud_T01.txt')
# diq_EB = np.loadtxt(f'../../Gosha/Barions/Goshe/Mu/newParametrs/04scalar_diq_ud00_EB.txt')
# diq_p1 = np.loadtxt(f'../../Gosha/Barions/Goshe/Mu/newParametrs/04scalar_diq_ud00.txt')
diq = np.loadtxt(f'../../Gosha/Barions/Goshe/04axial_diq_us35.txt')
# T_EB = diq_EB[:, 0]
# T_p1 = diq_p1[:, 0]
# T = diq_p1[:, 0]
# T_q = mdiq_ud_A[:, 0]
# T_q_2 = mdiq_ud[:, 0]
# m1_diq_A = mdiq_ud_A[:, 1]
# m2_diq_A = mdiq_ud_A[:, 2]
# m1_diq_B = mdiq_ud_B[:, 1]
# m2_diq_B = mdiq_ud_B[:, 2]
# m1_diq = mdiq_ud[:, 1]
# m2_diq = mdiq_ud[:, 2]

# # m_diq_EB = diq_EB[:, 1]
# # m_diq_p1 = diq_p1[:, 1]
# # m_diq = diq[:, 1]
# u_q = mueq_u_eb[:, 2]
# s_q = mueq_s_eb[:, 2]
# d_q = mueq_u_eb[:, 2]
#
# d_q_p1 = mueq_u_p1[:, 2]
# u_q_p1 = mueq_d_p1[:, 2]
# s_q_p1 = mueq_s_1[:, 2]



# u_q_m = u_q[:, 2]
# d_q_m = d_q[:, 2]
# s_q_m = s_q[:, 2]

# ax.plot(T_EB, m_diq_EB)
# ax.plot(T_p1, m_diq_p1)
# ax.plot(T, m_diq)
# ax.plot(T_q, u_q_eb, label="$A(m_u)$", linewidth=3, color="blue")
# ax.plot(T_q, s_q_eb, '--', label="$A(m_s)$", linewidth=3, color="blue")
# ax.plot(T_q, d_q_p1, label="$B(m_u)$", linewidth=3, color="orange")
# ax.plot(T_q, s_q_p1, '--', label="$B(m_s)$", linewidth=3, color="orange" )
ax.plot(diq[:, 0], diq[:, 1], label="$m_1A$", linewidth=1, color="green")
ax.plot(diq[:, 0],diq[:, 1], label="$m_1B$", linewidth=1, color="blue")
# ax.plot(u_q[:, 0], u_q_m, label="$m_1A$", linewidth=1, color="green")
# ax.plot(s_q[:, 0], s_q_m, label="$m_1B$", linewidth=1, color="blue")
# ax.plot(T_q, m2_diq_B, '--', label="$m_2B$", linewidth=1,  color="blue")
# ax.plot(T_q_2, m1_diq, label="$m_1$", linewidth=1, color="black")
# ax.plot(T_q_2, m2_diq, '--', label="$m_2$", linewidth=1,  color="black")
ax.set(xlim=[0.003, 0.450])
ax.tick_params(which='major', direction='in', width=2, length=6)
for axis in ['top', 'bottom', 'left', 'right']:
    ax.spines[axis].set_linewidth(2)
plt.legend()
plt.xlabel("$T$, GeV", fontsize=18, )
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.ylabel("masses, GeV", fontsize=18)
plt.show()
# fig.savefig('mass_q.eps')

# for i in mu:
#     lambda_bar = np.loadtxt(f'../../Gosha/Barions/Goshe/04scalar_diq_ud{i}.txt') if int(i) < 26 else np.loadtxt(f'../../../Gosha/Barions/Goshe/Mu/')
#     print(i)
#     T_l = lambda_bar[:, 0]
#     l_b = lambda_bar[:, 1]
#     Td_List.append(l_b[0])
#     #      mq_u2_mq_s = lambda_bar[:, 2]
#     # ax.plot(T_l, l_b, label=f'proton {i}')
# ax.plot(mu, Td_List, label=f'porog {i}')
# ax.legend()
# plt.show()

# pd_Gconst = np.loadtxt('./pd_Gconst.txt')
#
# for i in mu:
#     proton = np.loadtxt(f'../../Gosha/Barions/Goshe/04proton{i}.txt')
#     lambda_bar = np.loadtxt(f'../../Gosha/Barions/Goshe/04lambda_bar{i}.txt')
#     psi = np.loadtxt(f'../../Gosha/Barions/Goshe/04psi_bar{i}.txt')
#
#     T_p = proton[:, 0]
#     M_p = proton[:, 1]
#     mu3 = proton[:, 2]
#     Td_List.append(T_p[np.argmin(np.absolute(M_p - mu3))])
#     print(Td_List)
#     T_l = lambda_bar[:, 0]
#     l_b = lambda_bar[:, 1]
#     mq_u2_mq_s = lambda_bar[:, 2]
#     T_lambda_bar_list.append(T_l[np.argmin(np.absolute(l_b - mq_u2_mq_s))])
#
#     T_psi = psi[:, 0]
#     psi_bar = psi[:, 1]
#     mq_u2_mq_s_psi = psi[:, 2]
#     T_psi_list.append(T_l[np.argmin(np.absolute(psi_bar - mq_u2_mq_s_psi))])
#
# print(Td_List)
# print(T_lambda_bar_list)
#
# mu_2 = [0.00, 0.01, 0.03, 0.05, 0.07, 0.09, 0.11, 0.13, 0.15, 0.17, 0.19, 0.21, 0.23, 0.25, 0.26, 0.27, 0.28, 0.29,
#         0.30, 0.31, 0.32, 0.33, 0.34, 0.35]
# np.savetxt('T_mu_lambda.txt', np.c_[mu_2, T_lambda_bar_list], delimiter='    ')
# np.savetxt('T_mu_protons.txt', np.c_[mu_2, Td_List], delimiter='    ')
# np.savetxt('T_mu_psi.txt', np.c_[mu_2, T_psi_list], delimiter='    ')
#
# plt.scatter(mu_2, Td_List, c="g", label="proton")
# plt.scatter(mu_2, T_lambda_bar_list, c="r", label="lambda", marker="^")
# plt.scatter(mu_2, T_psi_list, c="b", label="psi", marker="*")
# plt.plot(pd_Gconst[:, 0], pd_Gconst[:, 1], label='pd_Gconst', linewidth=3)
#
# plt.xlabel("mu")
# plt.ylabel("T(MeV)")
# plt.legend()
# plt.show()
