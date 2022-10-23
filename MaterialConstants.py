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
        m0 = 9.1 * 10 ** (-31)  # electron mass in dormancy [kg]

        '''
        THE PARAMETERS OF QUANTUM WELL ARE DESCRIBED BELLOW:
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            The heterostructure has two QWs 6 and 3 nm wide.
            QW material - GaAs
        '''
        self.WIDTH_FIRST_QW = 6  # [nm]
        self.WIDTH_SECOND_QW = 3  # [nm]
        # LINK [http://www.matprop.ru/GaAs_bandstr]
        self.m_GaAs = 0.63 * m0  # effective electron mass in GaAs [kg]
        self.Eg_GaAs = 1.519 - 5.405 * 10 ** (-4) * (self.t ** 2) / (self.t + 204)  # energy gap GaAs at T K [eV]

        '''
        THE PARAMETERS OF BARRIERS ARE DESCRIBED BELLOW:
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            The barriers material is Al_(x)Ga_(1-x)As there x equal 30 %,  at that Al_(0.4)Ga_(0.6)As is between
        QW.
        '''
        self.MOLE_FRACTION = np.array([0.3, 0.4])  # x in Al_(x)Ga_(1-x)As
        self.WIDTH_BARRIERS = 2  # nm
        # Energy gap of barrier AlGaAs
        self.Eg_AlGaAs = np.array(1.424 + 1.247 * self.MOLE_FRACTION)  # link [http://www.matprop.ru/AlGaAs_basic]
        # effective electron mass in AlGaAs AT 300 K ! [kg]
        self.m_AlGaAs = np.array(0.063 + 0.083 * self.MOLE_FRACTION) * m0  # link [http://www.matprop.ru/AlGaAs_basic]
