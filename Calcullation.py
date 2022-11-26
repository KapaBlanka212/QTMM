import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import newton
from scipy.integrate import simpson, quad, fixed_quad

from MaterialConstants import StructureCONST
from Physics import Physics

CONST = StructureCONST(300)


def plot_t11_energy(steps, potential):
    U = potential
    min_E = U[0] / 100
    max_E = U[0] * 0.999
    root_e_find_space = np.geomspace(min_E, max_E, steps)
    t11_element_value = [Physics.find_energy(energy).real for energy in root_e_find_space]
    figure_e = plt.figure()
    ax = figure_e.add_subplot()
    ax.set_ylabel('Value element T[2,2] ', fontsize=12)
    ax.set_xlabel('E, eV', fontsize=12)
    ax.plot(root_e_find_space, t11_element_value)
    ax.plot(root_e_find_space, np.zeros(root_e_find_space.shape))
    plt.show()


if __name__ == "__main__":

    initial_guess = [0.07 * CONST.e, 0.16 * CONST.e, 0.28 * CONST.e]
    eqn_energy = lambda e: Physics.find_energy(e).real
    print(eqn_energy(0.0701 * CONST.e), eqn_energy(0.16 * CONST.e), eqn_energy(0.28 * CONST.e))
    E = []
    try:
        for x0 in initial_guess:
            energy_ans = newton(eqn_energy, x0)
            E.append(energy_ans)
    except:
        Exception
    print(E)
    steps = 1000
    z_min = -10 * 10 ** -9
    z_max = 15 * 10 ** -9
    z_space = np.linspace(z_min, z_max, steps)

    wave_functon = np.copy([Physics.find_wave_function(0.16 * CONST.e, item) for item in z_space])
    int = np.sqrt(quad(lambda x: np.abs(Physics.find_wave_function(0.16 * CONST.e, x)[0, :]) ** 2, z_min, z_max))
    print(int[0])

    figure_w = plt.figure()
    ax = figure_w.add_subplot()
    ax.set_ylabel('Wave Function', fontsize=12)
    ax.set_xlabel('z, nm', fontsize=12)
    ax.plot(z_space * 10 ** 9, wave_functon[:, 0], label='Full WF')
    ax.legend()
    plt.show()


