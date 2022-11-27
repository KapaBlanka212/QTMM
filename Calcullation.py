import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import newton
from scipy.integrate import simpson, quad, fixed_quad

from MaterialConstants import StructureCONST
from Physics import Physics

CONST = StructureCONST(300)


def plot_t11_energy(step, potential):
    U = potential
    min_E = U[0] / 100
    max_E = U[0] * 0.999
    root_e_find_space = np.geomspace(min_E, max_E, step)
    t11_element_value = [Physics.find_energy(energy)[1, 1].real for energy in root_e_find_space]
    figure_e = plt.figure()
    ax = figure_e.add_subplot()
    ax.set_ylabel('Value element T[2,2] ', fontsize=12)
    ax.set_xlabel('E, eV', fontsize=12)
    ax.plot(root_e_find_space, t11_element_value)
    ax.plot(root_e_find_space, np.zeros(root_e_find_space.shape))
    plt.show()


def plot_conductive_band(energy, potential):
    U0, U1 = potential
    x = np.array([0.0,
                  CONST.WIDTH_FIRST_QW,
                  CONST.WIDTH_BARRIERS + CONST.WIDTH_FIRST_QW,
                  CONST.WIDTH_SECOND_QW + CONST.WIDTH_FIRST_QW + CONST.WIDTH_BARRIERS,
                  CONST.WIDTH_SECOND_QW + CONST.WIDTH_FIRST_QW + CONST.WIDTH_BARRIERS])
    y = np.array([U0, 0.0, U1,   0.0, U0])

    fig = plt.figure()
    ax = fig.add_subplot()
    ax.step(x * 10 ** 9, y, linewidth=2, where='pre')
    for ans in energy:
        ax.plot([0.0 * 10 ** 9, (CONST.WIDTH_SECOND_QW + CONST.WIDTH_FIRST_QW + CONST.WIDTH_BARRIERS) * 10 ** 9],
                [ans, ans])
    ax.set_ylabel('Eg, eV', fontsize=15)
    ax.set_xlabel('x, nm', fontsize=15)
    plt.show()


def plot_wave_function(energy, z_array):
    wave_functions = [np.copy([Physics.find_wave_function(energy_i, z) for z in z_space]) for energy_i in energy]
    figure_e = plt.figure()
    ax = figure_e.add_subplot()
    ax.set_ylabel('Wave Function', fontsize=12)
    ax.set_xlabel('z, eV', fontsize=12)
    for wave_function_i in range(len(wave_functions)):
        ax.plot(z_array * 10 ** 9, wave_functions[wave_function_i][:, 0],
                label=f"Wave Function for E = {energy[wave_function_i]}")
    ax.plot(z_array * 10 ** 9, np.zeros(z_space.shape))
    ax.legend()
    plt.show()


def plot_norm_wave_function(energy, z_array):
    wave_functions_array = [np.copy([Physics.find_wave_function(energy_i, z) for z in z_space])
                            for energy_i in energy]
    wave_functions_sqr = [np.copy([np.abs(Physics.find_wave_function(energy_i, z)) ** 2 for z in z_space])
                          for energy_i in energy]
    integrals = [np.sqrt(simpson(wave_functions_sqr[i][:, 0], z_array, axis=0))
                 for i in range(len(wave_functions_sqr))]
    figure_e = plt.figure()
    ax = figure_e.add_subplot()
    ax.set_ylabel('Wave Function', fontsize=12)
    ax.set_xlabel('z, eV', fontsize=12)
    for wave_function_i in range(len(wave_functions_array)):
        ax.plot(z_array * 10 ** 9, wave_functions_array[wave_function_i][:, 0] / integrals[wave_function_i],
                label=f"Wave Function for E = {energy[wave_function_i]}")
    ax.plot(z_array * 10 ** 9, np.zeros(z_space.shape))
    ax.legend()
    plt.savefig("WF")
    plt.show()


def norm_wave_function(wave_function_not_norm, z):
    ints = np.sqrt(quad(wave_function_not_norm, z[0], z[-1]))
    return ints


if __name__ == "__main__":
    initial_guess = [0.06, 0.16, 0.28]
    steps = 1000
    z_space = np.linspace(-20, 20, steps) * 10 ** -9
    E = []
    steps = 1000
    for x0 in initial_guess:
        ans = newton(lambda x: Physics.find_energy(x)[1, 1].real, x0)
        E.append(ans)
    
    з
    plot_wave_function(E, z_space)
    plot_norm_wave_function(E, z_space)


