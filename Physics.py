import numpy as np
from cmath import exp
import matplotlib as pls


class Physics:
    def __init__(self):
        self.conduction_band = None

    @staticmethod
    def exp_body(k, a):
        out_body = 1j * k * a
        return out_body

    @staticmethod
    def qw_matrix(k, a):
        exp_body = Physics.exp_body(k, a)
        out_matrix = np.array([[exp(exp_body), 0],
                               [0, exp(-exp_body)]])
        return out_matrix

    @staticmethod
    def barrier_matrix(mass_, k):
        out_matrix = np.array([[1, 1],
                               [1j * k / mass_, 1j * k / mass_]])
        return out_matrix

    @staticmethod
    def inverse_barrier_matrix(mass_, k):
        out_matrix = np.array([[1 / 2, mass_ / (2j * k)],
                               [1 / 2, mass_ / (2j * k)]])
        return out_matrix

    def plot_conduction_band(self):

        pass
