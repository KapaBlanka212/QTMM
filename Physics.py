import numpy as np
from cmath import exp, sqrt

from MaterialConstants import StructureCONST
CONST = StructureCONST(300)


class Physics:

    @staticmethod
    def k_(mass_: float, energy: float, potential: float) -> float:
        energy *= CONST.e
        potential *= CONST.e
        k_out = sqrt(2 * mass_ * (energy - potential)) / CONST.h_
        return k_out

    @staticmethod
    def exp_body(k: float, a: float) -> float:
        width = a * 10 ** (-9)
        out_body = 1j * k * width
        return out_body

    @staticmethod
    def qw_matrix(a: float, k: float) -> np.ndarray:
        exp_body = Physics.exp_body(k, a)
        out_matrix = np.array([[exp(exp_body), 0.0],
                               [0.0, exp(-exp_body)]], dtype=object)
        return out_matrix

    @staticmethod
    def barrier_matrix(mass_: float, k: float) -> np.ndarray:
        out_matrix = np.array([[1, 1],
                               [(1j * k) / mass_, (-1j * k) / mass_]], dtype=object)
        return out_matrix

    @staticmethod
    def inverse_barrier_matrix(mass_: float, k: float) -> np.ndarray:
        out_matrix = 1 / 2 * np.array([[1, mass_ / (1j * k)],
                                       [1, - mass_ / (1j * k)]], dtype=object)
        return out_matrix

    @staticmethod
    def transfer_matrix(energy: float) -> float:
        # Describe barriers and K
        U0, U1 = np.array(CONST.Eg_GaAs - CONST.Eg_AlGaAs)
        k0 = k4 = Physics.k_(CONST.m_AlGaAs[0], energy, U0)
        k1 = k3 = Physics.k_(CONST.m_GaAs, energy, 0.0)
        k2 = Physics.k_(CONST.m_AlGaAs[1], energy, U1)

        m4_inv = Physics.inverse_barrier_matrix(CONST.m_AlGaAs[0], k4)
        m3 = Physics.barrier_matrix(CONST.m_GaAs, k3)
        n2 = Physics.qw_matrix(CONST.WIDTH_SECOND_QW, k3)
        m3_inv = Physics.inverse_barrier_matrix(CONST.m_GaAs, k3)
        m2 = Physics.barrier_matrix(CONST.m_AlGaAs[1], k2)
        n1 = Physics.qw_matrix(CONST.WIDTH_BARRIERS, k2)
        m2_inv = Physics.inverse_barrier_matrix(CONST.m_AlGaAs[1], k2)
        m1 = Physics.barrier_matrix(CONST.m_GaAs, k1)
        n0 = Physics.qw_matrix(CONST.WIDTH_FIRST_QW, k1)
        m1_inv = Physics.inverse_barrier_matrix(CONST.m_GaAs, k1)
        m0 = Physics.barrier_matrix(CONST.m_AlGaAs[0], k0)

        t = m4_inv @ m3 @ n2 @ m3_inv @ m2 @ n1 @ m2_inv @ m1 @ n0 @ m1_inv @ m0
        return t[1, 1]

