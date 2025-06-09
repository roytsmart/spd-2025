import pathlib
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.visualization
import named_arrays as na
import optika
import ccd_snr
from ._path import default_path


def eqe() -> pathlib.Path:

    ccd = ccd_snr.ccd()

    wavelength = na.geomspace(1, 10000, axis="wavelength", num=1001) << u.AA

    energy = wavelength.to(u.eV, equivalencies=u.spectral())

    eqe = ccd.quantum_efficiency_effective(
        rays=optika.rays.RayVectorArray(
            wavelength=wavelength,
            direction=na.Cartesian3dVectorArray(0, 0, 1),
        ),
        normal=na.Cartesian3dVectorArray(0, 0, -1),
    )

    # Store the wavelengths at which the QE was measured
    wavelength_measured = ccd.quantum_efficiency_measured.inputs

    # Store the QE measurements
    qe_measured = ccd.quantum_efficiency_measured.outputs

    with astropy.visualization.quantity_support():
        fig, ax = plt.subplots(
            figsize=(7, 6),
            constrained_layout=True,
        )
        ax2 = ax.twiny()
        ax2.invert_xaxis()
        color = "tab:red"
        na.plt.plot(wavelength, eqe, ax=ax, color=color);
        na.plt.plot(energy, eqe, ax=ax2, color=color);
        na.plt.scatter(wavelength_measured, qe_measured, ax=ax, color=color);
        ax.set_xscale("log");
        ax2.set_xscale("log");
        ax.set_xlabel(f"wavelength ({ax.get_xlabel()})");
        ax2.set_xlabel(f"energy ({ax2.get_xlabel()})", labelpad=16)
        ax.set_ylabel(f"effective quantum efficiency");

    path = default_path / "eqe.svg"

    fig.savefig(path)

    return path
