import numpy as np
from scipy.optimize import newton
import matplotlib.pyplot as plt

from MaterialConstants import StructureCONST
from Physics import Physics

CONST = StructureCONST(300)


if __name__ == "__main__":
    U0 = CONST.Eg_GaAs - CONST.Eg_AlGaAs[0]
    steps = 1000000
    min_E = 0.04
    max_E = 0.98 * U0
    energy_space = np.geomspace(min_E, max_E, steps)
    matrix_element_value = [[(Physics.transfer_matrix(item)).real, Physics.transfer_matrix(item).imag]
                            for item in energy_space]
    fig1 = plt.figure()
    ax1 = fig1.add_subplot()
    ax1.set_ylabel('Value element T[2,2] ', fontsize=15)
    ax1.set_xlabel('E, eV', fontsize=15)
    ax1.plot(energy_space, matrix_element_value)
    plt.show()
    '''
    PLOT UNIT:
    ~~~~~~~~~~
    '''
    """
    potential = np.array([CONST.Eg_GaAs - CONST.Eg_AlGaAs[0], CONST.Eg_GaAs - CONST.Eg_AlGaAs[1]])
    x = np.array([0.0,
                  CONST.WIDTH_FIRST_QW,
                  CONST.WIDTH_BARRIERS + CONST.WIDTH_FIRST_QW,
                  CONST.WIDTH_SECOND_QW + CONST.WIDTH_FIRST_QW + CONST.WIDTH_BARRIERS,
                  CONST.WIDTH_SECOND_QW + CONST.WIDTH_FIRST_QW + CONST.WIDTH_BARRIERS])
    y = np.array([CONST.Eg_GaAs - CONST.Eg_AlGaAs[0],
                  0.0,
                  CONST.Eg_GaAs - CONST.Eg_AlGaAs[1],
                  0.0,
                  CONST.Eg_GaAs - CONST.Eg_AlGaAs[0]])

    fig = plt.figure()
    ax = fig.add_subplot()
    ax.step(x, y, linewidth=2, where='pre')
    ax.set_ylabel('Eg, eV', fontsize=15)
    ax.set_xlabel('x, nm', fontsize=15)
    ax.set(xlim=[0, np.max(x) + 1], xticks=np.linspace(np.min(x), np.max(x) + 1, 13),
           ylim=[np.min(y), np.max(y) + 0.5], yticks=np.linspace(np.min(y), np.max(y), 5))

    plt.savefig("Eg.png")
    """

