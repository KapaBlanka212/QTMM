import numpy as np


class StructureCONST:
    """
    This class keep constant of material that create heterostructure GaAs / Al_(x)Ga_(1 - x)As
    """
    def __init__(self, temperature: float = 300.0):
        """
        :param temperature: K
        """
        super(StructureCONST, self).__init__()
        self.t = temperature
        m0 = 9.11 * 10 ** (-31)  # electron mass in dormancy [kg]
        self.h_ = 1.05 * 10 ** (-34)
        self.e = 1.6 * 10 ** (-19)

        '''
        THE PARAMETERS OF QUANTUM WELL ARE DESCRIBED BELLOW:
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            The heterostructure has two QWs 6 and 3 nm wide.
            QW material - GaAs
        '''

        self.WIDTH_FIRST_QW = 6 * 10 ** -9  # [m]
        self.WIDTH_SECOND_QW = 3 * 10 ** -9  # [m]
        # LINK [http://www.matprop.ru/GaAs_bandstr]

        self.m_GaAs = 0.063 * m0  # effective electron mass in GaAs [kg]
        self.Eg_GaAs = 4.07 * self.e  # electron affinity at T K [eV]

        '''
        THE PARAMETERS OF BARRIERS ARE DESCRIBED BELLOW:
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            The barriers material is Al_(x)Ga_(1-x)As there x equal 30 %,  at that Al_(0.4)Ga_(0.6)As is between
        QW.
        '''
        self.MOLE_FRACTION = np.array([0.4, 0.3])  # x in Al_(x)Ga_(1-x)As
        self.WIDTH_BARRIERS = 2 * 10 ** -9  # m

        # Electron affinity of barrier AlGaAs  LINK [http://www.matprop.ru/AlGaAs_basic]
        self.Eg_AlGaAs = (4.07 - 1.1 * self.MOLE_FRACTION) * self.e

        # Effective electron mass in AlGaAs AT 300 K ! [kg]
        self.m_AlGaAs = (0.063 + 0.083 * self.MOLE_FRACTION) * m0  # LINK [http://www.matprop.ru/AlGaAs_basic]

        """
        CALCULATE POTENTIAL
        """
        self.U = self.Eg_GaAs - self.Eg_AlGaAs

