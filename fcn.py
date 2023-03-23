import numpy as np

from IntegA import integral_A
from const import G, Gv, a0, a1, T0, a2, a3, m_u0, Nc, pi, K, m_d0, b3, mu_u0, b4, mu_d0, mu_s0, m_s0
from forMu import intmu
from integaphi import  integaphi
from integphi import integphi


def fcn(x,  T):
    f= np.zeros(8)
    mu_u = x[5]
    mu_d = x[6]
    mu_s = x[7]

    GT = G * (1.0 - 0.2 * x[3] * x[4] - 0.2 * (x[3] ** 3 + x[4] ** 3))
    GvT = Gv * (1.0 - 0.2 * x[3] * x[4] - 0.2 * (x[3] ** 3 + x[4] ** 3))

    ia_u = integral_A(x[0], mu_u, x[3], x[4], T)
    ia_d = integral_A(x[1], mu_d, x[3], x[4], T)
    ia_s = integral_A(x[2], mu_s, x[3], x[4], T)

    b2 = a0 + a1 * (T0 / T) + a2 * (
            T0 / T) ** 2 + a3 * (T0 / T) ** 3

    f[0] = 1.0 - m_u0 / x[0] + GT * Nc / pi ** 2 * ia_u - \
           K * Nc ** 2 / 8.0 / pi ** 4 * x[1] * x[2] / x[
               0] * ia_d * ia_s

    f[1] = 1.0 - m_d0 / x[1] + GT * Nc / pi ** 2 * ia_d - \
           K * Nc ** 2 / 8.0 / pi ** 4 * x[0] * x[2] / x[
               1] * ia_u * ia_s

    f[2] = 1.0 - m_s0 / x[2] + GT * Nc / pi ** 2 * ia_s - \
           K * Nc ** 2 / 8.0 / pi ** 4 * x[0] * x[1] / x[
               2] * ia_u * ia_d

    f[3] = (-b2 * x[4] - b3 * x[3] * x[3] + b4 * x[3] * x[4] * x[4]) / 2.0 - \
           2.0 * Nc / 2.0 / pi ** 2 / T ** 3 * integphi(x[0], mu_u, x[3], x[4], T) - \
           2.0 * Nc / 2.0 / pi ** 2 / T ** 3 * integphi(x[1], mu_d, x[3], x[4], T) - \
           2.0 * Nc / 2.0 / pi ** 2 / T ** 3 * integphi(x[2], mu_s, x[3], x[4], T)

    f[4] = (-b2 * x[3] - b3 * x[4] * x[4] + b4 * x[4] * x[3] * x[3]) / 2.0 - \
           2.0 * Nc / 2.0 / pi ** 2 / T ** 3 * integaphi(x[0], mu_u, x[3], x[4], T) - \
           2.0 * Nc / 2.0 / pi ** 2 / T ** 3 * integaphi(x[1], mu_d, x[3], x[4], T) - \
           2.0 * Nc / 2.0 / pi ** 2 / T ** 3 * integaphi(x[2], mu_s, x[3], x[4], T)

    f[5] = x[5] - mu_u0 + 2.0 * Nc / pi ** 2 * GvT * (
            intmu(x[0], x[5], x[3], x[4], T) + intmu(x[1], x[6], x[3], x[4], T) + intmu(x[2], x[7], x[3], x[4], T))

    f[6] = x[6] - mu_d0 + 2.0 * Nc / pi ** 2 * GvT * (
            intmu(x[0], x[5], x[3], x[4], T) + intmu(x[1], x[6], x[3], x[4], T,) + intmu(x[2], x[7], x[3], x[4], T))

    f[7] = x[7] - mu_s0 +  2.0 * Nc / pi ** 2 * GvT * (
            intmu(x[0], x[5], x[3], x[4], T) + intmu(x[1], x[6], x[3], x[4], T) + intmu(x[2], x[7], x[3], x[4], T))
    return f
