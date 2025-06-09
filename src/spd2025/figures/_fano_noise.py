import pathlib
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.visualization
import named_arrays as na
import ccd_snr
from ._path import default_path


def fano_noise() -> pathlib.Path:

    ccd = ccd_snr.ccd()

    wavelength = na.geomspace(1, 10000, axis="wavelength", num=1001) << u.AA

    energy = wavelength.to(u.eV, equivalencies=u.spectral())

    fano_factor = ccd.fano_factor(wavelength)

    with astropy.visualization.quantity_support():
        fig, ax = plt.subplots(
            figsize=(7, 6),
            constrained_layout=True,
        )
        ax2 = ax.twiny()
        ax2.invert_xaxis()
        na.plt.plot(wavelength, fano_factor, ax=ax, color="mediumslateblue");
        na.plt.plot(energy, fano_factor, ax=ax2, color="mediumslateblue");
        ax.set_xscale("log");
        ax2.set_xscale("log");
        # ax.set_yscale("log");
        # ax2.set_yscale("log");
        ax.set_xlabel(f"wavelength ({ax.get_xlabel()})");
        ax2.set_xlabel(f"energy ({ax2.get_xlabel()})", labelpad=16)
        ax.set_ylabel(f"Fano factor ({ax.get_ylabel()})");

    path = default_path / "fano_factor.svg"

    fig.savefig(path)

    return path
