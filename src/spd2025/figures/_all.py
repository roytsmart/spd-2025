import pathlib
import matplotlib as mpl
from ._quantum_yield import quantum_yield
from ._schematic import schematic
from ._penetration_depth import penetration_depth
from ._differential_cce import differential_cce
from ._cce import cce
from ._absorbance import absorbance
from ._eqe import eqe
from ._fano_noise import fano_noise
from ._vmr_1_simple import vmr_simple
from ._vmr_2_shot import vmr_shot
from ._vmr_3_fano import vmr_fano
from ._vmr_4_pcc import vmr_pcc
from ._vmr_5_total import vmr_total
from ._vmr_ratio import vmr_ratio
from ._charge_diffusion import charge_diffusion

__all__ = [
    "all",
]


def all() -> list[pathlib.Path]:
    """
    Generate all the figures
    """

    factor = 2

    settings = {
        "font.size": 22,
        "font.family": "Calibri",
        "lines.linewidth": 1.5 * factor,
        "patch.linewidth": 1 * factor,
        "xtick.major.size": 3.5 * factor,
        "xtick.minor.size": 2 * factor,
        "xtick.major.width": 0.8 * factor,
        "xtick.minor.width": 0.6 * factor,
        "xtick.major.pad": 3.5 * factor,
        "xtick.minor.pad": 3.4 * factor,
        "ytick.major.size": 3.5 * factor,
        "ytick.minor.size": 2 * factor,
        "ytick.major.width": 0.8 * factor,
        "ytick.minor.width": 0.6 * factor,
        "ytick.major.pad": 3.5 * factor,
        "ytick.minor.pad": 3.4 * factor,
    }

    with mpl.rc_context(settings):
        paths = [
            quantum_yield(),
            schematic(),
            differential_cce(),
            penetration_depth(),
            cce(),
            absorbance(),
            eqe(),
            fano_noise(),
            vmr_simple(),
            vmr_shot(),
            vmr_fano(),
            vmr_pcc(),
            vmr_total(),
            vmr_ratio(),
            charge_diffusion(),
        ]

    return paths


if __name__ == '__main__':

    import spd2025

    spd2025.figures.all()
