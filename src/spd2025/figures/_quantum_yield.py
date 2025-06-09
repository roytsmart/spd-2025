import pathlib
import matplotlib.pyplot as plt
import astropy.units as u
import named_arrays as na
import optika
from ._path import default_path


def quantum_yield() -> pathlib.Path:

    # Define an array of wavelengths
    wavelength = na.geomspace(1, 10000, axis="wavelength", num=1001) << u.AA

    energy = wavelength.to(u.eV, equivalencies=u.spectral())

    # Compute the quantum yield
    iqy = optika.sensors.quantum_yield_ideal(wavelength)

    # Plot the quantum yield vs wavelength
    fig, ax = plt.subplots(
        figsize=(7, 6),
        constrained_layout=True,
    )
    ax2 = ax.twiny()
    ax2.invert_xaxis()
    na.plt.plot(wavelength, iqy, ax=ax);
    na.plt.plot(energy, iqy, ax=ax2);
    ax.set_xscale("log");
    ax2.set_xscale("log");
    ax.set_yscale("log");
    ax2.set_yscale("log");
    ax.set_xlabel(f"wavelength ({wavelength.unit:latex_inline})");
    ax2.set_xlabel(f"energy ({energy.unit:latex_inline})", labelpad=16)
    ax.set_ylabel(f"quantum yield ({iqy.unit:latex_inline})");

    path = default_path / "quantum_yield.svg"

    fig.savefig(path)

    return path
