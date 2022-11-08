from typing import List, Union, Tuple, NoReturn
import numpy as np
from cmath import exp, sqrt
import matplotlib.pyplot as plt

from MaterialConstants import StructureCONST
CONST = StructureCONST(300)


class Physics:

    @staticmethod
    def k_in_qw(m: float, energy: float):
        energy *= CONST.e
        k_out = sqrt(2 * m * energy) / CONST.h_
        return k_out

    @staticmethod
    def k_in_barrier(m: float, energy: float, potential: float):
        energy *= CONST.e
        potential *= CONST.e
        k_out = sqrt(2 * m * (energy - potential)) / CONST.h_
        return k_out

    @staticmethod
    def exp_body(k, a):
        a *= 10 ** -9
        out_body = 1j * k * a
        return out_body

    @staticmethod
    def qw_matrix(mass_, energy: float, a, k: float):
        exp_body = Physics.exp_body(k, a)
        out_matrix = np.array([[exp(exp_body), 0],
                               [0, exp(-exp_body)]])
        return out_matrix

    @staticmethod
    def barrier_matrix(mass_, k):
        out_matrix = np.array([[1, 1],
                               [1j * k / mass_, -1j * k / mass_]])
        return out_matrix

    @staticmethod
    def inverse_barrier_matrix(mass_, k):
        out_matrix = np.array([[1 / 2, - 1j * mass_ / (2 * k)],
                               [1 / 2, 1j * mass_ / (2 * k)]])
        return out_matrix

    @staticmethod
    def transfer_matrix(energy: float):
        # Describe barriers and wave vectors
        U0, U1 = np.array(CONST.Eg_GaAs - CONST.Eg_AlGaAs)
        k0 = k4 = Physics.k_in_barrier(CONST.m_AlGaAs[0], energy, U0)
        k1 = k3 = Physics.k_in_qw(CONST.m_GaAs, energy)
        k2 = Physics.k_in_barrier(CONST.m_AlGaAs[1], energy, U1)

        # Describe matrices
        m4_inv = Physics.inverse_barrier_matrix(CONST.m_AlGaAs[0], k4)
        m3 = Physics.barrier_matrix(CONST.m_GaAs, k3)
        n2 = Physics.qw_matrix(CONST.m_GaAs, energy,
                               CONST.WIDTH_FIRST_QW + CONST.WIDTH_SECOND_QW + CONST.WIDTH_BARRIERS, k1)
        m3_inv = Physics.inverse_barrier_matrix(CONST.m_GaAs, k3)
        m2 = Physics.barrier_matrix(CONST.m_AlGaAs[1], k2)
        n3 = Physics.qw_matrix(CONST.m_AlGaAs[1], energy, CONST.WIDTH_FIRST_QW + CONST.WIDTH_BARRIERS, k2)
        m2_inv = Physics.inverse_barrier_matrix(CONST.m_AlGaAs[1], k2)
        m1 = Physics.barrier_matrix(CONST.m_GaAs, k1)
        n1 = Physics.qw_matrix(CONST.m_GaAs, energy, CONST.WIDTH_FIRST_QW, k1)
        m1_inv = Physics.inverse_barrier_matrix(CONST.m_GaAs, k1)
        m0 = Physics.barrier_matrix(CONST.m_AlGaAs[0], k0)

        t = m4_inv @ m3 @ n2 @ m3_inv @ m2 @ n3 @ m2_inv @ m1 @ n1 @ m1_inv @ m0
        return t[1, 1]

