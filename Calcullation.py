import numpy as np
from scipy.optimize import newton
import matplotlib.pyplot as plt

from MaterialConstants import StructureCONST
from Physics import Physics

CONST = StructureCONST(300)


if __name__ == "__main__":
    U0, U1 = np.array(CONST.Eg_GaAs - CONST.Eg_AlGaAs)
    steps = 10000
    min_E = U0 / 1000
    max_E = 0.99 * U0

    energy_space = np.hstack([np.geomspace(min_E, max_E / 2, steps), np.linspace(max_E / 2, max_E, steps)])
    wave_f = [Physics.find_wave_function(x, 0.07) for x in x_space]
    matrix_element_value = [np.real(Physics.transfer_matrix(item)) for item in energy_space]

    fig1 = plt.figure()
    ax1 = fig1.add_subplot()
    ax1.set_ylabel('Value element T[2,2] ', fontsize=15)
    ax1.set_xlabel('E, eV', fontsize=15)
    ax1.plot(energy_space, matrix_element_value)
    plt.show()
    initial_guees = [0.04, 0.16, 0.29]

    ans_list = []
    for x0 in initial_guees:
        ans = newton(Physics.transfer_matrix, x0)
        ans_list.append(ans.real)
    print(ans_list)
    '''
    PLOT UNIT:
    ~~~~~~~~~~
    '''

    x = np.array([0.0,
                  CONST.WIDTH_FIRST_QW,
                  CONST.WIDTH_BARRIERS + CONST.WIDTH_FIRST_QW,
                  CONST.WIDTH_SECOND_QW + CONST.WIDTH_FIRST_QW + CONST.WIDTH_BARRIERS,
                  CONST.WIDTH_SECOND_QW + CONST.WIDTH_FIRST_QW + CONST.WIDTH_BARRIERS])
    y = np.array([U0,
                  0.0,
                  U1,
                  0.0,
                  U0])
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.step(x, y, linewidth=2, where='pre')
    ax.set_ylabel('Eg, eV', fontsize=15)
    ax.set_xlabel('x, nm', fontsize=15)
    ax.set(xlim=[-1, np.max(x) + 1], xticks=np.linspace(np.min(x), np.max(x) + 1, 13),
           ylim=[np.min(y), np.max(y) + 0.1], yticks=np.linspace(np.min(y), np.max(y), 10))
    for answer in ans_list:
        ax.plot([x[0], x[-1]], [answer, answer])
    plt.show()
