import numpy as np
import pandas as pd
import scipy

from scipy import integrate

mu_u, mu_d, mu_s = 0.0, 0.0, 0.0
mq_u, mq_d, mq_s = 0.0, 0.0, 0.0
T = 0.0
mq1, mu_q1, mq2, mu_q2 = 0.0, 0.0, 0.0, 0.0
lam, E01, E02 = 0.0, 0.0, 0.0
phi, aphi = 0.0, 0.0
mq, mu_q = 0.0, 0.0
Nc = 3.0
Nf = 3.0
pi = 3.14
G = 4.3
Gv = 0.0 * G
Gdiq = 0.705 * G
K = 89.8
L = 0.652

x = np.zeros(2)
x[0] = 1.1
x[1] = 0.1

mueq_u_mass_phi = pd.read_csv('mueq05Gv0GT_u_mass_phi_mu0.txt', header=None)
mueq_dmass = pd.read_csv('mueq05Gv0GT_dmass_mu0.txt', header=None)
mueq_smass = pd.read_csv('mueq05Gv0GT_smass_mu0.txt', header=None)

mq_u_list = mueq_u_mass_phi[2]
phi_list = mueq_u_mass_phi[3]
aphi_list = mueq_u_mass_phi[4]
mq_d_list = mueq_dmass[2]
mu_d_list = mueq_dmass[3]

T_list = mueq_smass[0]
mq_s_list = mueq_smass[2]
mu_s_list = mueq_smass[3]


def fermi1(en, mu, phi, aphi):
    fer1 = np.exp(-(en - mu) / T)
    u1 = 1.0 + 3.0 * (phi + aphi * fer1) * fer1 + fer1 ** 3

    return ((phi + 2.0 * aphi * fer1) * fer1 + fer1 ** 3) / u1


def fermi2(en, mu, phi, aphi):
    fer2 = np.exp(-(en + mu) / T)
    u2 = 1.0 + 3.0 * (aphi + phi * fer2) * fer2 + fer2 ** 3
    return ((aphi + 2.0 * phi * fer2) * fer2 + fer2 ** 3) / u2


def HeavisideTheta(arg):
    if (arg < 0.0):
        return 0.0
    else:
        return 1.0


def fB4_1(x):
    f1 = fermi2(x, mu_q2, phi, aphi)
    p = x ** 2 - mq2 ** 2
    if p < 0.0:
        p = 0.0
    if f1< 0.0001:
        f1=0
    return np.sqrt(p) * (1.0 - f1) / 2.0 / lam / (x - E02)



def fB3_1(x):
    f1 = fermi2(x, mu_q2, phi, aphi)
    p = x ** 2 - mq2 ** 2
    if p < 0.0:
        p = 0.0
    return np.sqrt(p) * f1 / 2.0 / lam / (x + E02)



def fB2_1(x):
    f1 = fermi2(x, mu_q1, phi, aphi)
    p = x ** 2 - mq1 ** 2
    if p < 0.0:
        p = 0.0
    return np.sqrt(p) * (1.0 - f1) / 2.0 / lam / (x + E01)




def fB1_1(x):
    f1 = fermi1(x, mu_q1, phi, aphi)
    p = x ** 2 - mq1 ** 2
    if p < 0.0:
        p = 0.0
    return np.sqrt(p) * f1 / 2.0 / lam / (x - E01)




def int_ReB1(mmq1, mmu_q1, mmq2, mmu_q2):
    global mq1, mu_q1, mq2, mu_q2
    mq1 = mmq1
    mu_q1 = mmu_q1
    mq2 = mmq2
    mu_q2 = mmu_q2

    if mu_q2 ** 2 - mq2 ** 2 + mq1 ** 2 > 0.0:
        a = np.sqrt(
            mu_q2 ** 2 - mq2 ** 2 + mq1 ** 2)
    else:
        a = mq1

    b = np.sqrt(L ** 2 + mq1 ** 2)

    if E01 < 0.0:
        return 4.0 * integrate.quad(fB1_1, a, b)[0]
    else:
        c = E01
        if c == a:
            return 4.0 * integrate.quad(fB1_1, a, b, points=[c], )[0]
        else:
            return 4.0 * integrate.quad(fB1_1, a, b, points=[c], )[0]


def int_ReB2(mmq1, mmu_q1, mmq2, mmu_q2):
    global mq1, mu_q1, mq2, mu_q2
    mq1 = mmq1
    mu_q1 = mmu_q1
    mq2 = mmq2
    mu_q2 = mmu_q2

    if mu_q2 ** 2 - mq2 ** 2 + mq1 ** 2 > 0.0:
        a = np.sqrt(
            mu_q2 ** 2 - mq2 ** 2 + mq1 ** 2)
    else:
        a = mq1

    b = np.sqrt(L ** 2 + mq1 ** 2)

    if E01 > 0.0:
        return 4.0 * integrate.quad(fB2_1, a, b)[0]
    else:
        c = -E01
        if c == a:
            return 4.0 * integrate.quad(fB2_1, a, b, points=[a, b])[0]
        else:
            return 4.0 * integrate.quad(fB2_1, a, b, points=[c])[0]


def int_ReB3(mmq1, mmu_q1, mmq2, mmu_q2):
    global mq1, mu_q1, mq2, mu_q2
    mq1 = mmq1
    mu_q1 = mmu_q1
    mq2 = mmq2
    mu_q2 = mmu_q2

    a = mu_q2
    b = np.sqrt(L ** 2 + mu_q2 ** 2)

    if E02 > 0.0:
        return 4.0 * integrate.quad(fB3_1, a, b)[0]
    else:
        c = -E02
        if c == a:
            return 4.0 * integrate.quad(fB3_1, a, b, points=[a, b])[0]
        else:
            return 4.0 * integrate.quad(fB3_1, a, b, points=[c])[0]


def int_ReB4(mmq1, mmu_q1, mmq2, mmu_q2):
    global mq1, mu_q1, mq2, mu_q2
    mq1 = mmq1
    mu_q1 = mmu_q1
    mq2 = mmq2
    mu_q2 = mmu_q2

    a = mu_q2
    b = np.sqrt(L ** 2 + mu_q2 ** 2)

    if E02 < 0.0:
        return 4.0 * integrate.quad(fB4_1, a, b)[0]
    else:
        c = E02
        if c == a:
            return 4.0 * integrate.quad(fB4_1, a, b, points=[a, b])[0]
        else:
            return 4.0 * integrate.quad(fB4_1, a, b, points=[c], epsrel=0.001)[0]


def fmq(x):
    en = np.sqrt(mq * mq + x * x)
    return x ** 2 / en * (1.0 - fermi1(en, mu_q, phi,
                                       aphi) - fermi2(en,
                                                      mu_q,
                                                      phi,
                                                      aphi))


def integral_A(mmq, mmu_q):
    global mq, mu_q
    mq = mmq
    mu_q = mmu_q

    if mu_u ** 2 - mq_u ** 2 > 0.0:
        a = np.sqrt(mu_u ** 2 - mq_u ** 2)
    else:
        a = 0.0
    return -4.0 * integrate.quad(fmq, a, L)[0]


def integral_ImB(mmes, mmq1, mmu_q1, mmq2, mmu_q2):
    global mq1, mu_q1, mq2, mu_q2, lam, E01, E02
    Im_B1, Im_B2, Im_B3, Im_B4 = 0.0, 0.0, 0.0, 0.0
    mes = mmes
    mq1 = mmq1
    mu_q1 = mmu_q1
    mq2 = mmq2
    mu_q2 = mmu_q2

    lam = mes - mu_q2 + mu_q1
    E01 = (lam ** 2 + mq1 ** 2 - mq2 ** 2) / 2.0 / lam
    E02 = (lam ** 2 - mq1 ** 2 + mq2 ** 2) / 2.0 / lam

    Theta = HeavisideTheta(E01 - mq1) * HeavisideTheta(
        np.sqrt(L ** 2 + mq1 ** 2) - E01)

    if Theta > 0.0:
        Im_B1 = 4.0 * pi * np.sqrt(E01 ** 2 - mq1 ** 2) * fermi1(
            E01,
            mu_q1,
            phi,
            aphi) * Theta / 2.0 / lam

    Theta = HeavisideTheta(-E01 - mq1) * HeavisideTheta(
        np.sqrt(L ** 2 + mq1 ** 2) + E01)

    if Theta > 0.0:
        Im_B2 = 4.0 * pi * np.sqrt(E01 ** 2 - mq1 ** 2) * fermi1(
            E01,
            mu_q1,
            phi,
            aphi) * Theta / 2.0 / lam

    Theta = HeavisideTheta(-E02 - mq2) * HeavisideTheta(
        np.sqrt(L ** 2 + mq2 ** 2) + E02)

    if Theta > 0.0:
        Im_B3 = 4.0 * pi * np.sqrt(E02 ** 2 - mq2 ** 2) * (
                1 - fermi2(E02,
                           mu_q2,
                           phi,
                           aphi)) * Theta / 2.0 / lam

    Theta = HeavisideTheta(E02 - mq2) * HeavisideTheta(
        np.sqrt(L ** 2 + mq2 ** 2) - E02)

    if Theta > 0.0:
        Im_B4 = 4.0 * pi * np.sqrt(E02 ** 2 - mq2 ** 2) * (
                1 - fermi2(E02,
                           mu_q2,
                           phi,
                           aphi)) * Theta / 2.0 / lam
    # print(f'Im_B1 + Im_B2 - Im_B3 - Im_B4= {Im_B1 + Im_B2 - Im_B3 - Im_B4}')
    return Im_B1 + Im_B2 - Im_B3 - Im_B4


def integral_ReB(mmes, mmq1, mmu_q1, mmq2, mmu_q2):
    global mq1, mu_q1, mq2, mu_q2, lam, E01, E02
    mes = mmes
    mq1 = mmq1
    mu_q1 = mmu_q1
    mq2 = mmq2
    mu_q2 = mmu_q2

    lam = mes - mu_q2 + mu_q1
    E01 = (lam ** 2 + mq1 ** 2 - mq2 ** 2) / 2.0 / lam
    E02 = (lam ** 2 - mq1 ** 2 + mq2 ** 2) / 2.0 / lam

    Re_B1 = int_ReB1(mq1, mu_q1, mq2,
                     mu_q2)

    Re_B2 = int_ReB2(mq1, mu_q1, mq2,
                     mu_q2)

    Re_B3 = int_ReB3(mq1, mu_q1, mq2,
                     mu_q2)

    Re_B4 = int_ReB4(mq1, mu_q1, mq2,
                     mu_q2)

    print(Re_B1, Re_B2, Re_B3, Re_B4)

    return -Re_B1 - Re_B2 + Re_B3 + Re_B4


def fcn(x):
    f = np.zeros(2)
    # print(f'x-{x}')
    ia_u = integral_A(mq_u, mu_u)
    # ia_d = integral_A(mq_d, mu_d)
    # ia_s = integral_A(mq_s, mu_s)
    # ia_s_min = integral_A(mq_s, - mu_s)
    ia_d_min = integral_A(mq_d, -mu_d)
    Gdiq_f = Gdiq / 4.0

    B0_ud = integral_ReB(x[0], mq_u, mu_u, mq_d,
                         -mu_d)
    # B0_us = integral_ReB(x[0], mq_u, mu_u, mq_s, -mu_s)

    Im_B0ud = integral_ImB(x[0], mq_u, mu_u, mq_d,
                           -mu_d)
    # Im_B0_us = integral_ImB(x[0], mq_u, mu_u, mq_s,
    #                        -mu_s)
    # print(f'Im_B0ud-{Im_B0ud}')
    # print(f'B0_ud-{B0_ud}')
    modul_du = B0_ud ** 2 + Im_B0ud ** 2
    # print(f'modul_du-{modul_du}')
    #
    # Pol_ds = - 2.0 / pi ** 2 * (ia_u + ia_s_min
    #                             + (mq_u ** 2 + mq_s ** 2 - 4.0 * mq_u * mq_s
    #                                - (x[0] + mu_u + mu_s) ** 2) * B0_us)

    f[0] = Gdiq_f * (x[0] + mu_u + mu_d) * x[1]
    + (1 - Gdiq_f * (ia_u + ia_d_min)) * Im_B0ud / modul_du

    f[1] = Gdiq_f * ((mu_u ** 2 + mu_d ** 2 - 4.0 * mu_u * mu_d
                      - (x[0] + mu_u + mu_d) ** 2) + 1.0 / 4.0 * x[1] ** 2)
    + (1 - Gdiq_f * (ia_u + ia_d_min)) * B0_ud / modul_du

    return f


for i in range(350):
    mq_u = mq_u_list[i]
    phi = phi_list[i]
    aphi = aphi_list[i]
    mq_d = mq_d_list[i]
    mu_d = mu_d_list[i]
    mq_s = mq_s_list[i]
    mu_s = mu_s_list[i]
    T = T_list[i]
    xguess = x.copy()
    x = scipy.optimize.root(fcn, xguess)['x']
    print(T, x[0])
