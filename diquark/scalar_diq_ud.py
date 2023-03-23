import math

import numpy as np
import pandas as pd
import scipy
from scipy import integrate


def HeavisideTheta(arg):
    if (arg < 0.0):
        return 0.0
    else:
        return 1.0


class Diquark:
    def __init__(self, mu_u, mu_d, mu_s, mq_u, mq_d, mq_s, phi, aphi, T):
        self.mu_u = mu_u
        self.mu_d = mu_d
        self.mu_s = mu_s
        self.mq_u = mq_u
        self.mq_d = mq_d
        self.mq_s = mq_s
        self.mq = 0.0
        self.mq1 = 0.0
        self.mq2 = 0.0
        self.mu_q = 0.0
        self.mu_q1 = 0.0
        self.mu_q2 = 0.0
        self.lam = 0.0
        self.E01 = 0.0
        self.E02 = 0.0
        self.phi = phi
        self.aphi = aphi
        self.T = T
        self.Nc = 3.0
        self.Nf = 3.0
        self.pi = 3.14
        self.G = 4.3
        self.Gv = 0.0
        self.Gdiq = 2.9706
        self.K = 89.8
        self.L = 0.652

    def fermi1(self, en, mu, phi, aphi):
        fer1 = math.exp(-(en - mu) / self.T)
        u1 = 1.0 + 3.0 * (phi + aphi * fer1) * fer1 + fer1 ** 3

        return ((phi + 2.0 * aphi * fer1) * fer1 + fer1 ** 3) / u1

    def fermi2(self, en, mu, phi, aphi):
#        print(self.T, en, mu, phi, aphi)
        fer2 = math.exp(-(en + mu) / self.T)
        u2 = 1.0 + 3.0 * (aphi + phi * fer2) * fer2 + fer2 ** 3

        return ((aphi + 2.0 * phi * fer2) * fer2 + fer2 ** 3) / u2

    def fcn(self, x):
        f = np.zeros(2)
        ia_u = self.integral_A(self.mq_u, self.mu_u)
        # ia_d = self.integral_A(self.mq_d, self.mu_d)
        # ia_s = self.integral_A(self.mq_s, self.mu_s)
        # ia_s_min = self.integral_A(self.mq_d, -self.mu_d)
        ia_d_min = self.integral_A(self.mq_s, -self.mu_s)

        Gdiq_f = self.Gdiq

        # print(self.mq_u, self.mu_u, self.mq_d, -self.mu_d)


        B0_ud = self.integral_Reb(x[0], self.mq_u, self.mu_u, self.mq_d, -self.mu_d)
        ImB0_ud = self.integral_ImB(x[0], self.mq_u, self.mu_u, self.mq_d, -self.mu_d)

        modul_ud = B0_ud ** 2 + ImB0_ud ** 2

        print(B0_ud, ImB0_ud)

        fmass_ud = ((self.mq_u - self.mq_d) ** 2 - (x[0] + self.mu_u + self.mu_d) ** 2
                    + x[1] ** 2 / 4.0)

        f[0] = 2.0 * Gdiq_f / self.pi ** 2 * (x[0] + self.mu_u + self.mu_d) * x[1]

        + (1.0 + 2.0 * Gdiq_f / self.pi ** 2 * (ia_u + ia_d_min)) * ImB0_ud / modul_ud

        f[1] = 2.0 * Gdiq_f / self.pi ** 2 * fmass_ud

        + (1.0 + 2.0 * Gdiq_f / self.pi ** 2 * (ia_u + ia_d_min)) * B0_ud / modul_ud

        return f

    def integral_A(self, mq, mu_q):
        self.mq = mq
        self.mu_q = mu_q
        p = self.mu_u ** 2 - self.mq_u ** 2
        if p > 0:
            a = math.sqrt(p)
        else:
            a = 0.0

        def fmq(x):
            en = math.sqrt(self.mq ** 2 + x ** 2)
            return x ** 2 / en * (
                    1.0 - self.fermi1(en, self.mu_q, self.phi, self.aphi) - self.fermi2(en, self.mu_q, self.phi,
                                                                                        self.aphi))

        return -4.0 * integrate.quad(fmq, a, self.L)[0]

    def integral_ImB(self, mdiq, mmq1, mmu_q1, mmq2, mmu_q2):

        self.mq1 = mmq1
        self.mu_q1 = mmu_q1
        self.mq2 = mmq2
        self.mu_q2 = mmu_q2
        self.lam = mdiq - self.mu_q2 + self.mu_q1

#

        self.E01 = (self.lam ** 2 + self.mq1 ** 2 - self.mq2 ** 2) / 2.0 / self.lam
        self.E02 = (self.lam ** 2 - self.mq1 ** 2 + self.mq2 ** 2) / 2.0 / self.lam

        # print(self.mq1, self.mu_q1, self.mq2, self.mu_q2, mdiq)
        Theta = HeavisideTheta((self.E01 - self.mq1) * (math.sqrt(self.L ** 2 + self.mq1 ** 2) - self.E01))
        Im_B1 = 0.0
        Im_B2 = 0.0
        Im_B3 = 0.0
        Im_B4 = 0.0
        if Theta > 0:
            Im_B1 = 4.0 * self.pi * math.sqrt(self.E01 ** 2 - self.mq1 ** 2) * self.fermi1(self.E01, self.mu_q1,
                                                                                           self.phi,
                                                                                           self.aphi)
        Theta = HeavisideTheta(-self.E01 - self.mq1) * HeavisideTheta(math.sqrt(self.L ** 2 + self.mq1 ** 2) + self.E01)

        if Theta > 0.0:
            Im_B2 = 4.0 * self.pi * math.sqrt(self.E01 ** 2 - self.mq1 ** 2) * self.fermi1(self.E01, self.mu_q1,
                                                                                           self.phi,
                                                                                           self.aphi) * Theta / 2.0 / self.lam

        Theta = HeavisideTheta(-self.E02 - self.mq2) * HeavisideTheta(
            math.sqrt(self.L ** 2 + self.mq2 ** 2) + self.E02)

        if Theta > 0.0:
            Im_B3 = 4.0 * self.pi * math.sqrt(self.E02 ** 2 - self.mq2 ** 2) * (
                    1.0 - self.fermi2(self.E02, self.mu_q2, self.phi, self.aphi)) * Theta / 2.0 / self.lam

        Theta = HeavisideTheta(self.E02 - self.mq2) * HeavisideTheta(math.sqrt(self.L ** 2 + self.mq2 ** 2) - self.E02)

        if Theta > 0.0:
            Im_B4 = 4.0 * self.pi * math.sqrt(self.E02 ** 2 - self.mq2 ** 2) * (
                    1.0 - self.fermi2(self.E02, self.mu_q2, self.phi, self.aphi)) * Theta / 2.0 / self.lam

  #      print(mdiq, Im_B1, Im_B2, Im_B3, Im_B4)
        return - Im_B1 - Im_B2 + Im_B3 + Im_B4

    def integral_Reb(self, mdiq, mq1, mu_q1, mq2, mu_q2):
        self.mq1 = mq1
        self.mu_q1 = mu_q1
        self.mq2 = mq2
        self.mu_q2 = mu_q2
        self.lam = mdiq - self.mu_q2 + self.mu_q1
        self.E01 = (self.lam ** 2 + self.mq1 ** 2 - self.mq2 ** 2) / 2.0 / self.lam
        self.E02 = (self.lam ** 2 - self.mq1 ** 2 + self.mq2 ** 2) / 2.0 / self.lam

        # print(self.mq1, self.mu_q1, self.mq2, self.mu_q2, mdiq)

        Re_B1 = self.int_ReB1()
        Re_B2 = self.int_ReB2()
        Re_B3 = self.int_ReB3()
        Re_B4 = self.int_ReB4()

#        print(Re_B1, Re_B2, Re_B3, Re_B4)
        return - Re_B1 - Re_B2 + Re_B3 + Re_B4

    def int_ReB1(self):
        a = self.mq1
        b = math.sqrt(self.L ** 2 + self.mq1 ** 2)

        def fb1(x):
            f1 = self.fermi2(x, self.mu_q1, self.phi, self.aphi)
            p = x ** 2 - self.mq1 ** 2
            if p < 0.0:
                p = 0.0
            return math.sqrt(p) * f1 / 2.0 / self.lam / (x - self.E01)

        if self.E01 < 0.0:
            return 4.0 * integrate.quad(fb1, a, b)[0]
        else:
            c = self.E01
            if c == a:
                return 4.0 * integrate.quad(fb1, a, b, points=[a, b])[0]
            else:
                return 4.0 * integrate.quad(fb1, a, b, points=[c])[0]

    def int_ReB2(self):
        a = self.mq1
        b = math.sqrt(self.L ** 2 + self.mq1 ** 2)

        def fb2(x):
            f1 = self.fermi2(x, self.mu_q1, self.phi, self.aphi)
            p = x ** 2 - self.mq1 ** 2
            if p < 0.0:
                p = 0.0
 #           print(self.lam)

            return math.sqrt(p) * (1.0 - f1) / 2.0 / self.lam / (x + self.E01)

        if self.E01 > 0.0:
            return 4.0 * integrate.quad(fb2, a, b)[0]
        else:
            c = -self.E01
            if c == a:
                return 4.0 * integrate.quad(fb2, a, b, points=[a, b])[0]
            else:
                return 4.0 * integrate.quad(fb2, a, b, points=[c])[0]

    def int_ReB3(self):
        a = self.mq2
        b = math.sqrt(self.L ** 2 + self.mq2 ** 2)

        def fb3(x):
            f1 = self.fermi2(x, self.mu_q2, self.phi, self.aphi)
            p = x ** 2 - self.mq2 ** 2
            if p < 0.0:
                p = 0.0
            return math.sqrt(p) * f1 / 2.0 / self.lam / (x + self.E02)

        if self.E02 < 0.0:
            return 4.0 * integrate.quad(fb3, a, b)[0]
        else:
            c = -self.E02
            if c == a:
                return 4.0 * integrate.quad(fb3, a, b, points=[a, b])[0]
            else:
                return 4.0 * integrate.quad(fb3, a, b, points=[c])[0]

    def int_ReB4(self):
        a = self.mq2
        b = math.sqrt(self.L ** 2 + self.mq2 ** 2)
        def fb4(x):

            f1 = self.fermi2(x, self.mu_q2, self.phi, self.aphi)
#            print(f1, self.mq2)
            p = x ** 2 - self.mq2 ** 2

            if p < 0.0:
                p = 0.0
            return math.sqrt(p) * (1.0 - f1) / 2.0 / self.lam
        print(self.E02)
        if self.E02 < 0.0:
            return 4.0 * integrate.quad(fb4, a, b)[0]
        else:
            c = self.E02
            if c == a:
                return 4.0 * integrate.quad(fb4, a, b, points=[a, b])[0]
            else:
                return 4.0 * integrate.quad(fb4, a, b, points=[c])[0]


x = np.zeros(2)
x[0] = 1.1
x[1] = 0.000001

mueq_u_mass_phi = np.loadtxt('./04mueq_u_mass_phi_mu00.txt')
mueq_dmass = np.loadtxt('./04mueq_dmass_mu00.txt')
mueq_smass = np.loadtxt('./04mueq_smass_mu00.txt')

mq_u_arr = mueq_u_mass_phi[:, 2]
phi_arr = mueq_u_mass_phi[:, 3]
aphi_arr = mueq_u_mass_phi[:, 4]

mq_d_arr = mueq_dmass[:, 2]
mu_d_arr = mueq_dmass[:, 3]

T_arr = mueq_smass[:, 0]
mq_s_arr = mueq_smass[:, 2]
mu_s_arr = mueq_smass[:, 3]

for i in range(350):
    diq = Diquark(mu_d_arr[i],
                  mu_d_arr[i],
                  mu_s_arr[i],
                  mq_u_arr[i],
                  mq_d_arr[i],
                  mq_s_arr[i],
                  phi_arr[i],
                  aphi_arr[i],
                  T_arr[i],
                  )
    xguess = x.copy()
    x = scipy.optimize.root(diq.fcn, xguess)['x']
    print(diq.T, x[0], x[1])

