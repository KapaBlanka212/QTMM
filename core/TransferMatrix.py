from dataclasses import dataclass

import numpy as np
from cmath import exp, sqrt
from typing import List

from core.sovler_setting import Setting
from MaterialConstants import StructureCONST
CONST = StructureCONST()


class WaveVector:
    def __init__(self, setting: Setting):
        self.setting = setting

    @property
    def wave_vector(self):
        return self.__calculation_wave_vector()

    def __calculation_wave_vector(self) -> float:
        k_out = sqrt(2 * self.setting.mass * (self.setting.energy - self.setting.potential)) / CONST.h_
        return k_out


class PropagationMatrix:
    def __init__(self, setting: Setting):
        # settings
        if not isinstance(setting, Setting):
            raise AttributeError('Setting not instance')
        self.setting = setting
        self.wave_vector = WaveVector(setting)

    @property
    def propagation_matrix(self):
        if self.setting.offset is not None:            return self.__calculate_offset_propagation_matrix()
        else:
            if self.setting.well_width is None:
                raise AttributeError('Well width must be not None')
            return  self.__calculate_propagation_matrix()

    def __calculate_propagation_matrix(self):
        out_matrix = np.array([
            [exp(1j * self.wave_vector.wave_vector * self.setting.well_width), 0.0],
            [0.0, exp(-1j * self.wave_vector.wave_vector * self.setting.well_width)]], np.dtype(np.complex128))
        return out_matrix

    def __calculate_offset_propagation_matrix(self):
        out_matrix = np.array([
            [exp(1j * self.wave_vector.wave_vector * self.setting.offset), 0.0],
            [0.0, exp(-1j * self.wave_vector.wave_vector * self.setting.offset)]], np.dtype(np.complex128))
        return out_matrix

class InterfaceMatrix:
    def __init__(self, setting: Setting):
        self.setting = setting
        self.wave_vector = WaveVector(self.setting)

    @property
    def inv_interface_matrix(self):
        return self.__calculate_inverse_interface_matrix()

    @property
    def interface_matrix(self):
        return self.__calculate_interface_matrix()

    def __calculate_interface_matrix(self) -> np.ndarray:
        elements = 1j * self.wave_vector.wave_vector / self.setting.mass
        out_matrix = np.array([
            [1, 1],
            [elements, -elements]], np.dtype(np.complex128))
        return out_matrix

    def __calculate_inverse_interface_matrix(self) -> np.ndarray:
        elements = 1j * self.setting.mass / self.wave_vector.wave_vector
        out_matrix = np.array([
            [1, -elements],
            [1, elements]], np.dtype(np.complex128))
        return 1 / 2 * out_matrix


class TransferMatrix:
    def __init__(self, layer: Setting):
        self.propagation_matrix = PropagationMatrix(layer)
        self.interface_matrix = InterfaceMatrix(layer)
        self.__transfer_matrix = None

    @property
    def transfer_matrix(self):
        self.__transfer_matrix = self.__calculate_transfer_matrix()
        return self.__transfer_matrix

    def __calculate_transfer_matrix(self):
        propagation_matrix = self.propagation_matrix.propagation_matrix
        inv_interface_matrix = self.interface_matrix.inv_interface_matrix
        interface_matrix = self.interface_matrix.interface_matrix
        matrix = interface_matrix @ propagation_matrix @ inv_interface_matrix
        return matrix

class FullTransferMatrix:
    def __init__(self, structure: List[Setting]):
        self.structure = structure

    @property
    def matrix(self):
        return self.__calculate_full_transfer_matrix()

    @property
    def matrix_elements(self):
        return self.matrix[1, 1]

    def __calculate_full_transfer_matrix(self):
        full_matrix = np.eye(2)
        for i in reversed(range(1, len(self.structure) - 1)):
            matrix = TransferMatrix(self.structure[i])
            full_matrix = full_matrix @ matrix.transfer_matrix
        inv_interface_matrix = InterfaceMatrix(self.structure[-1]).inv_interface_matrix
        interface_matrix = InterfaceMatrix(self.structure[0]).interface_matrix
        full_matrix = inv_interface_matrix @ full_matrix @ interface_matrix
        return full_matrix



