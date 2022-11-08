import numpy as np
from scipy.optimize import newton
import matplotlib.pyplot as plt

from MaterialConstants import StructureCONST
from Physics import Physics

CONST = StructureCONST(300)


if __name__ == "__main__":
    initial_guess = np.linspace(0, CONST.Eg_GaAs - CONST.Eg_AlGaAs[0], 100)
    ans_list = []
    for x0 in initial_guess:
        try:
            ans = newton(Physics.transfer_matrix, x0)
            if abs(np.imag(ans)) < 10 ** -20:  # Energy must be real :3
                if abs(Physics.transfer_matrix(ans.real)) < 10 ** -5:  # Check how close the answer is to 0
                    ans_list.append(np.round(ans.real, 5))
        except Exception:
            pass
    ans_array = np.copy(ans_list)

    '''
    PLOT UNIT:
    ~~~~~~~~~~
    '''
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
    ax.plot([np.min(x), np.max(x)], [ans_array, ans_array], linewidth=1)
    plt.show()
    plt.savefig("Eg")
    print(f'Energy level is {np.sort(np.transpose(np.unique(ans_array)))}')

