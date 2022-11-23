import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import newton
from scipy.integrate import simpson

from MaterialConstants import StructureCONST
from Physics import Physics

CONST = StructureCONST(300)


def plot_t11_energy(steps, potential):
    U = potential
    min_E = U[1] / 100
    max_E = U[1] * 0.999
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
    initial_guess = [0.04, 0.15, 0.28]
    eqn_energy = lambda e: Physics.find_energy(e).real
    E = []
    for x0 in initial_guess:
        energy_ans = newton(eqn_energy, x0)
        E.append(energy_ans)

    steps = 100
    z_min = -1 * 10 ** -9
    z_max = 8.99 * 10 ** -9
    z_space = np.linspace(z_min, z_max, steps)
    wave_functon = np.copy(Physics.find_wave_function(E[0], z_space))
    wave_functon1 = np.copy(Physics.find_wave_function(E[1], z_space))
    wave_functon2 = np.copy(Physics.find_wave_function(E[2], z_space))

    figure_f = plt.figure()
    ax = figure_f.add_subplot()
    ax.set_ylabel('Diff Wave Function', fontsize=12)
    ax.set_xlabel('z, nm', fontsize=12)
    ax.plot(z_space * 10 ** 9, wave_functon[:, 1])
    ax.plot(z_space * 10 ** 9, wave_functon1[:, 1])
    ax.plot(z_space * 10 ** 9, wave_functon2[:, 1])
    ax.plot(z_space * 10 ** 9, np.zeros(z_space.shape))

    figure_w = plt.figure()
    ax = figure_w.add_subplot()
    ax.set_ylabel('Wave Function', fontsize=12)
    ax.set_xlabel('z, nm', fontsize=12)
    ax.plot(z_space * 10 ** 9, wave_functon[:, 0])
    ax.plot(z_space * 10 ** 9, wave_functon1[:, 0])
    ax.plot(z_space * 10 ** 9, wave_functon2[:, 0])
    ax.plot(z_space * 10 ** 9, np.zeros(z_space.shape))
    plt.show()


