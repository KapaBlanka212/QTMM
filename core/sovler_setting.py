import warnings
from dataclasses import dataclass, fields
import numpy as np

from MaterialConstants import StructureCONST
CONST = StructureCONST()


@dataclass()
class Setting:
    _energy: np.typing.NDArray | float = None
    _potential: np.typing.NDArray | float = None
    _effective_mass: np.typing.NDArray | float = None
    _well_width: np.typing.NDArray | float = None
    _z0: np.typing.NDArray | float = None
    _z: np.typing.NDArray | float = None

    @classmethod
    def from_dict(cls, data: dict):
        field_mapping = {
            "mass": "_effective_mass",
            "well_width": "_well_width",
            "energy": "_energy",
            "potential": "_potential"
        }

        kwargs = {}
        for key, value in data.items():
            field_name = field_mapping.get(key, f"_{key}")
            if not any(f.name == field_name for f in fields(cls)):
                warnings.warn(f"Игнорируем неизвестный ключ: {key}", UserWarning)
                continue

            kwargs[field_name] = value

        return cls(**kwargs)

    @property
    def mass(self):
        return self._effective_mass

    @property
    def energy(self):
        return self._energy * CONST.e

    @property
    def potential(self):
        return self._potential * CONST.e

    @property
    def offset(self):
        if self._z is not None and self._z0 is not None:
            return self._z - self._z0
        return None

    @property
    def well_width(self):
        return self._well_width