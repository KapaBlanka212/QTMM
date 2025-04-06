import numpy as np
import scipy as sp

from core.TransferMatrix import FullTransferMatrix, Setting, CONST


class Solver:
    def __init__(self, structure_parameters: list[dict],
                 initial_guess: np.typing.NDArray | float):
        self.structure_parameters = structure_parameters
        self.initial_guess = initial_guess

    @property
    def answer(self):
        return np.real(np.array(self.solve_energy_equation()))

    def create_solver_setting(self, energy):
        setting_list: list = []
        for ELEMENT in self.structure_parameters:
            ELEMENT['energy'] = energy
            setting_list.append(Setting.from_dict(ELEMENT))
        return setting_list

    def energy_equation(self, energy):
        setting_list: list[Setting] = self.create_solver_setting(energy)
        transfer_matrix = FullTransferMatrix(setting_list)
        return transfer_matrix.matrix_elements

    def solve_energy_equation(self):
        answer_list = []
        for guess in self.initial_guess:
            answer_list.append(sp.optimize.newton(self.energy_equation, guess, maxiter=200))
        return answer_list

    def wave_function_equation(self):
        pass





# old realisation

# @staticmethod
# def find_energy(energy: float) -> np.ndarray:
#     # there we find E - energy eigenvalues
#     U = MaterialConstants.StructureCONST.U  # U0 = 0.44 eV, U1 = 0.33 eV
#
#     """
#     DESCRIBED K
#     """
#     k0 = TransferMatrix(energy=energy, potential=U[0]).k_(MaterialConstants.StructureCONST.m_AlGaAs[0])
#     k1 = TransferMatrix(energy=energy, potential=0).k_(MaterialConstants.StructureCONST.m_GaAs)
#     k2 = TransferMatrix(energy=energy, potential=U[1]).k_(MaterialConstants.StructureCONST.m_AlGaAs[1])
#     k3 = TransferMatrix(energy=energy, potential=0).k_(MaterialConstants.StructureCONST.m_GaAs)
#     k4 = TransferMatrix(energy=energy, potential=U[0]).k_(MaterialConstants.StructureCONST.m_AlGaAs[0])
#
#     """
#     DESCRIBED TRANSFER MATRICES
#     """
#
#     m0 = TransferMatrix.m_matrix(MaterialConstants.StructureCONST.m_AlGaAs[0], k0)
#     m1 = TransferMatrix.m_matrix(MaterialConstants.StructureCONST.m_GaAs, k1)
#     m1_inv = TransferMatrix.inverse_m_matrix(MaterialConstants.StructureCONST.m_GaAs, k1)
#     m2 = TransferMatrix.m_matrix(MaterialConstants.StructureCONST.m_AlGaAs[1], k2)
#     m2_inv = TransferMatrix.inverse_m_matrix(MaterialConstants.StructureCONST.m_AlGaAs[1], k2)
#     m3 = TransferMatrix.m_matrix(MaterialConstants.StructureCONST.m_GaAs, k3)
#     m3_inv = TransferMatrix.inverse_m_matrix(MaterialConstants.StructureCONST.m_GaAs, k3)
#     m4_inv = TransferMatrix.inverse_m_matrix(MaterialConstants.StructureCONST.m_AlGaAs[0], k4)
#
#     n1 = TransferMatrix.n_matrix(MaterialConstants.StructureCONST.WIDTH_FIRST_QW, k1)
#     n2 = TransferMatrix.n_matrix(MaterialConstants.StructureCONST.WIDTH_BARRIERS, k2)
#     n3 = TransferMatrix.n_matrix(MaterialConstants.StructureCONST.WIDTH_SECOND_QW, k3)
#
#     """
#     CALCULATED ELEMENT T11 OF THE FULL TRANSFER MATRIX
#     """
#
#     t1 = n1 @ m1_inv @ m0
#     t2 = n2 @ m2_inv @ m1
#     t3 = n3 @ m3_inv @ m2
#     full_t = m4_inv @ m3 @ n3 @ m3_inv @ m2 @ n2 @ m2_inv @ m1 @ n1 @ m1_inv @ m0
#     return full_t
#
#
# @staticmethod
# def find_wave_function(energy: float, z: float) -> List:
#     U = MaterialConstants.StructureCONST.U  # U0 = 0.44 eV, U1 = 0.33 eV
#
#     z0 = 0.0
#     z1 = MaterialConstants.StructureCONST.WIDTH_FIRST_QW
#     z2 = z1 + MaterialConstants.StructureCONST.WIDTH_BARRIERS
#     z3 = z2 + MaterialConstants.StructureCONST.WIDTH_SECOND_QW
#
#     """
#     DESCRIBED K
#     ~~~~~~~~~~~
#     """
#     k0 = TransferMatrix(energy=energy, potential=U[0]).k_(MaterialConstants.StructureCONST.m_AlGaAs[0])
#     k1 = TransferMatrix(energy=energy, potential=0).k_(MaterialConstants.StructureCONST.m_GaAs)
#     k2 = TransferMatrix(energy=energy, potential=U[1]).k_(MaterialConstants.StructureCONST.m_AlGaAs[1])
#     k3 = TransferMatrix(energy=energy, potential=0).k_(MaterialConstants.StructureCONST.m_GaAs)
#     k4 = TransferMatrix(energy=energy, potential=U[0]).k_(MaterialConstants.StructureCONST.m_AlGaAs[0])
#
#     """
#     DESCRIBED TRANSFER MATRICES
#     ~~~~~~~~~~~~~~~~~~~~~~~~~~~
#     """
#     m0 = TransferMatrix.m_matrix(MaterialConstants.StructureCONST.m_AlGaAs[0], k0)
#     m1 = TransferMatrix.m_matrix(MaterialConstants.StructureCONST.m_GaAs, k1)
#     m1_inv = TransferMatrix.inverse_m_matrix(MaterialConstants.StructureCONST.m_GaAs, k1)
#     m2 = TransferMatrix.m_matrix(MaterialConstants.StructureCONST.m_AlGaAs[1], k2)
#     m2_inv = TransferMatrix.inverse_m_matrix(MaterialConstants.StructureCONST.m_AlGaAs[1], k2)
#     m3 = TransferMatrix.m_matrix(MaterialConstants.StructureCONST.m_GaAs, k3)
#     m3_inv = TransferMatrix.inverse_m_matrix(MaterialConstants.StructureCONST.m_GaAs, k3)
#     m4 = TransferMatrix.m_matrix(MaterialConstants.StructureCONST.m_AlGaAs[0], k4)
#     m4_inv = TransferMatrix.inverse_m_matrix(MaterialConstants.StructureCONST.m_AlGaAs[0], k4)
#
#     n1_1 = TransferMatrix.n_matrix(MaterialConstants.StructureCONST.WIDTH_FIRST_QW, k1)
#     n2_1 = TransferMatrix.n_matrix(MaterialConstants.StructureCONST.WIDTH_BARRIERS, k2)
#     n3_1 = TransferMatrix.n_matrix(MaterialConstants.StructureCONST.WIDTH_SECOND_QW, k3)
#
#     n1 = TransferMatrix.n_z_matrix(z, z0, k1)
#     n2 = TransferMatrix.n_z_matrix(z, z1, k2)
#     n3 = TransferMatrix.n_z_matrix(z, z2, k3)
#     n4 = TransferMatrix.n_z_matrix(z, z3, k4)
#
#     C0 = np.array([[0], [1]])
#     t1 = n1_1 @ m1_inv @ m0
#     t2 = n2_1 @ m2_inv @ m1
#     t3 = n3_1 @ m3_inv @ m2
#
#     C1 = n1 @ m1_inv @ m0 @ C0
#     C2 = n2 @ m2_inv @ m1 @ t1 @ C0
#     C3 = n3 @ m3_inv @ m2 @ t2 @ t1 @ C0
#     C4 = n4 @ m4_inv @ m3 @ t3 @ t2 @ t1 @ C0
#
#     if z < z0:
#         return m0 @ (C0 * exp(-1j * k0 * z))
#     elif z0 < z < z1:
#         return m1 @ C1
#     elif z1 < z < z2:
#         return m2 @ C2
#     elif z2 < z < z3:
#         return m3 @ C3
#     elif z > z3:
#         return m4 @ C4