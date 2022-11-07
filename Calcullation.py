import logging as log

import numpy as np
from scipy.optimize import newton

from MaterialConstants import StructureCONST
from Physics import Physics

CONST = StructureCONST(300)


if __name__ == "__main__":
    initual_guess = np.linspace(0, CONST.Eg_GaAs - CONST.Eg_AlGaAs[0], 100)
    ans_list = []
    for x0 in initual_guess:
        try:
            ans = newton(Physics.transfer_matrix, x0)
            if abs(np.imag(ans)) < 10 ** -19:
                ans_list.append(np.round(ans.real, 2))
        except Exception:
            pass
    ans_array = np.copy(ans_list)
    print(np.sort(np.transpose(np.unique(ans_array))))
