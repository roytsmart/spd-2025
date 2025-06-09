import pathlib
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.visualization
import named_arrays as na
import optika
import ccd_snr
from ._path import default_path


def absorbance() -> pathlib.Path:

    ccd = ccd_snr.ccd()

    wavelength = na.geomspace(1, 10000, axis="wavelength", num=1001) << u.AA

    energy = wavelength.to(u.eV, equivalencies=u.spectral())

    absorbance = ccd.absorbance(
        rays=optika.rays.RayVectorArray(
            wavelength=wavelength,
            direction=na.Cartesian3dVectorArray(0, 0, 1),
        ),
        normal=na.Cartesian3dVectorArray(0, 0, -1),
    ).average

    with astropy.visualization.quantity_support():
        fig, ax = plt.subplots(
            figsize=(7, 6),
            constrained_layout=True,
        )
        ax2 = ax.twiny()
        ax2.invert_xaxis()
        na.plt.plot(wavelength, absorbance, ax=ax, color="royalblue");
        na.plt.plot(energy, absorbance, ax=ax2, color="royalblue");
        ax.set_xscale("log");
        ax2.set_xscale("log");
        ax.set_xlabel(f"wavelength ({ax.get_xlabel()})");
        ax2.set_xlabel(f"energy ({ax2.get_xlabel()})", labelpad=16)
        ax.set_ylabel(f"absorbance");

    path = default_path / "absorbance.svg"

    fig.savefig(path)

    return path
