import numpy as np
from cmath import exp, sqrt
from typing import List

from MaterialConstants import StructureCONST
CONST = StructureCONST()


class Physics:

    @staticmethod
    def k_(mass_: float, energy: float, potential: float) -> float:
        k_out = sqrt(2 * mass_ * (energy - potential)) / CONST.h_
        return 1j * abs(k_out)

    @staticmethod
    def real_k(mass_: float, energy: float) -> float:
        k_out = np.sqrt(2 * mass_ * energy) / CONST.h_
        return k_out

    @staticmethod
    def exp_body(k: float, a: float) -> complex:
        out_body = 1j * k * a
        return out_body

    @staticmethod
    def n_matrix(a: float, k: float) -> np.ndarray:
        exp_body = Physics.exp_body(k, a)
        out_matrix = np.array([[exp(exp_body), 0.0],
                               [0.0, exp(-exp_body)]])
        return out_matrix

    @staticmethod
    def n_z_matrix(z: float, z0: float, k: float) -> np.ndarray:
        out_matrix = np.array([[exp(1j * k * (z - z0)), 0.0],
                               [0.0, exp(-1j * k * (z - z0))]])
        return out_matrix

    @staticmethod
    def m_matrix(mass_: float, k: float) -> np.ndarray:
        elements = 1j * k / mass_
        out_matrix = np.array([[1, 1],
                               [elements, -elements]])
        return out_matrix

    @staticmethod
    def inverse_m_matrix(mass_: float, k: float) -> np.ndarray:
        elements = 1j * mass_ / k
        out_matrix = 1 / 2 * np.array([[1, -elements],
                                       [1, elements]])
        return out_matrix

    @staticmethod
    def find_energy(energy: float) -> float:
        # there we find E - energy eigenvalues
        U = CONST.U  # U0 = 0.44 eV, U1 = 0.33 eV

        """
        DESCRIBED K
        """
        k0 = Physics.k_(CONST.m_AlGaAs[0], energy, U[0])
        k1 = Physics.real_k(CONST.m_GaAs, energy)
        k2 = Physics.k_(CONST.m_AlGaAs[1], energy, U[1])
        k3 = Physics.real_k(CONST.m_GaAs, energy)
        k4 = Physics.k_(CONST.m_AlGaAs[0], energy, U[0])

        """
        DESCRIBED TRANSFER MATRICES
        """

        m0 = Physics.m_matrix(CONST.m_AlGaAs[0], k0)
        m1 = Physics.m_matrix(CONST.m_GaAs, k1)
        m1_inv = Physics.inverse_m_matrix(CONST.m_GaAs, k1)
        m2 = Physics.m_matrix(CONST.m_AlGaAs[1], k2)
        m2_inv = Physics.inverse_m_matrix(CONST.m_AlGaAs[1], k2)
        m3 = Physics.m_matrix(CONST.m_GaAs, k3)
        m3_inv = Physics.inverse_m_matrix(CONST.m_GaAs, k3)
        m4_inv = Physics.inverse_m_matrix(CONST.m_AlGaAs[0], k4)

        n1 = Physics.n_matrix(CONST.WIDTH_FIRST_QW, k1)
        n2 = Physics.n_matrix(CONST.WIDTH_BARRIERS, k2)
        n3 = Physics.n_matrix(CONST.WIDTH_SECOND_QW, k3)

        """
        CALCULATED ELEMENT T11 OF THE FULL TRANSFER MATRIX 
        """

        t1 = n1 @ m1_inv @ m0
        t2 = n2 @ m2_inv @ m1
        t3 = n3 @ m3_inv @ m2
        full_t = m4_inv @ m3 @ t3 @ t2 @ t1

        return full_t[1, 1]

    @staticmethod
    def find_wave_function(energy: float, z: float) -> List:
        U = CONST.U  # U0 = 0.44 eV, U1 = 0.33 eV

        z0 = 0.0
        z1 = CONST.WIDTH_FIRST_QW
        z2 = z1 + CONST.WIDTH_BARRIERS
        z3 = z2 + CONST.WIDTH_SECOND_QW

        """
        DESCRIBED K
        ~~~~~~~~~~~
        """
        k0 = Physics.k_(CONST.m_AlGaAs[0], energy, U[0])
        k1 = Physics.real_k(CONST.m_GaAs, energy)
        k2 = Physics.k_(CONST.m_AlGaAs[1], energy, U[1])
        k3 = Physics.real_k(CONST.m_GaAs, energy)
        k4 = Physics.k_(CONST.m_AlGaAs[0], energy, U[0])

        """
        DESCRIBED TRANSFER MATRICES
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~
        """
        m0 = Physics.m_matrix(CONST.m_AlGaAs[0], k0)
        m1 = Physics.m_matrix(CONST.m_GaAs, k1)
        m1_inv = Physics.inverse_m_matrix(CONST.m_GaAs, k1)
        m2 = Physics.m_matrix(CONST.m_AlGaAs[1], k2)
        m2_inv = Physics.inverse_m_matrix(CONST.m_AlGaAs[1], k2)
        m3 = Physics.m_matrix(CONST.m_GaAs, k3)
        m3_inv = Physics.inverse_m_matrix(CONST.m_GaAs, k3)
        m4 = Physics.m_matrix(CONST.m_AlGaAs[0], k4)
        m4_inv = Physics.inverse_m_matrix(CONST.m_AlGaAs[0], k4)

        n1 = Physics.n_z_matrix(z, z0, k1)
        n2 = Physics.n_z_matrix(z, z1, k2)
        n3 = Physics.n_z_matrix(z, z2, k3)
        n4 = Physics.n_z_matrix(z, z3, k4)

        C0 = np.array([[0], [1]])
        C1 = n1 @ m1_inv @ m0 @ C0
        C2 = n2 @ m2_inv @ m1 @ C1
        C3 = n3 @ m3_inv @ m2 @ C2
        C4 = n4 @ m4_inv @ m3 @ C3
        C4[1] = 0.0
        if z < z0:
            return m0 @ (C0 * exp(-1j * k0 * z))
        elif z0 < z < z1:
            return m1 @ C1
        elif z1 < z < z2:
            return m2 @ C2
        elif z2 < z < z3:
            return m3 @ C3
        elif z > z3:
            return m4 @ C4


