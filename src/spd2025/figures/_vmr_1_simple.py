import pathlib
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.visualization
import named_arrays as na
import optika
import ccd_snr
from ._path import default_path


def vmr_simple() -> pathlib.Path:

    ccd = ccd_snr.ccd()

    wavelength = na.geomspace(1.5, 10000, axis="wavelength", num=1001) << u.AA

    energy = wavelength.to(u.eV, equivalencies=u.spectral())

    rays = optika.rays.RayVectorArray(
        wavelength=wavelength,
        direction=na.Cartesian3dVectorArray(0, 0, 1),
    )

    normal = na.Cartesian3dVectorArray(0, 0, -1)

    eqe = ccd.quantum_efficiency_effective(rays, normal)

    vmr_simple = 1 / eqe * u.photon

    color_simple = "gray"

    with astropy.visualization.quantity_support():
        fig, ax = plt.subplots(
            figsize=(7, 6),
            constrained_layout=True,
        )
        ax2 = ax.twiny()
        ax2.invert_xaxis()
        na.plt.plot(wavelength, vmr_simple, ax=ax, color=color_simple, label="simple");
        na.plt.plot(energy, vmr_simple, ax=ax2, color=color_simple);
        ax.set_xscale("log");
        ax2.set_xscale("log");
        ax.set_xlabel(f"wavelength ({ax.get_xlabel()})");
        ax2.set_xlabel(f"energy ({ax2.get_xlabel()})", labelpad=16)
        ax.set_ylabel(f"variance-to-mean ratio ({ax.get_ylabel()})");
        ax.set_ylim(bottom=-1, top=15)
        ax.legend()

    path = default_path / "vmr_simple.svg"

    fig.savefig(path)

    return path
