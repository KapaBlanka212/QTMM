import numpy as np

from MaterialConstants import StructureCONST
from core.solver import Solver
import matplotlib.pyplot as plt

CONST = StructureCONST()

def main():
    # tested structure AlGaAs/GaAs/AlGas/GaAs/AlGas
    potential = CONST.U
    initial_guess = np.eye(3).flatten() + 0.1

    algaas_0_1_dict = {'potential': potential[1],'mass': CONST.m_AlGaAs[1],
                       'well_width': CONST.WIDTH_BARRIERS}
    gaas_1 = {'potential': 0.0,'mass': CONST.m_GaAs,
              'well_width': CONST.WIDTH_FIRST_QW}

    structure_ = [algaas_0_1_dict, gaas_1, algaas_0_1_dict, gaas_1, algaas_0_1_dict, gaas_1, algaas_0_1_dict]
    structure = structure_
    solver = Solver(structure, initial_guess)
    plot_t11_energy(1000, potential, solver)
    #print(solver.answer)

def plot_t11_energy(step, potential, solver: Solver):
    U = potential
    min_E = U[0] / 100
    max_E = U[0] * 0.999
    root_e_find_space = np.geomspace(min_E, max_E, step)
    t11_element_value = [solver.energy_equation(energy).real for energy in root_e_find_space]
    figure_e = plt.figure()
    ax = figure_e.add_subplot()
    ax.set_ylabel('Value element T[2,2] ', fontsize=12)
    ax.set_xlabel('E, eV', fontsize=12)
    ax.plot(root_e_find_space, t11_element_value)
    ax.plot(root_e_find_space, np.zeros(root_e_find_space.shape))
    plt.show()



if __name__ == "__main__":
    main()
