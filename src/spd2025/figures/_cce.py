import pathlib
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.visualization
import named_arrays as na
import optika
import ccd_snr
from ._path import default_path


def cce() -> pathlib.Path:

    ccd = ccd_snr.ccd()

    wavelength = na.geomspace(1, 10000, axis="wavelength", num=1001) << u.AA

    energy = wavelength.to(u.eV, equivalencies=u.spectral())

    cce = ccd.charge_collection_efficiency(
        rays=optika.rays.RayVectorArray(
            wavelength=wavelength,
            direction=na.Cartesian3dVectorArray(0, 0, 1),
        ),
        normal=na.Cartesian3dVectorArray(0, 0, -1),
    )

    with astropy.visualization.quantity_support():
        fig, ax = plt.subplots(
            figsize=(7, 6),
            constrained_layout=True,
        )
        ax2 = ax.twiny()
        ax2.invert_xaxis()
        na.plt.plot(wavelength, cce, ax=ax, color="tab:purple");
        na.plt.plot(energy, cce, ax=ax2, color="tab:purple");
        ax.set_xscale("log");
        ax2.set_xscale("log");
        ax.set_xlabel(f"wavelength ({ax.get_xlabel()})");
        ax2.set_xlabel(f"energy ({ax2.get_xlabel()})", labelpad=16)
        ax.set_ylabel(f"charge-collection efficiency");

    path = default_path / "cce.svg"

    fig.savefig(path)

    return path
