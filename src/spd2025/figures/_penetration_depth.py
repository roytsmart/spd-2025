import pathlib
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.visualization
import named_arrays as na
import optika
import ccd_snr
from ._path import default_path


def penetration_depth() -> pathlib.Path:

    ccd = ccd_snr.ccd()

    wavelength = na.geomspace(1, 10000, axis="wavelength", num=1001) << u.AA

    energy = wavelength.to(u.eV, equivalencies=u.spectral())

    depth = 1 / optika.chemicals.Chemical("Si").absorption(wavelength)

    depth = depth.to(u.nm)

    with astropy.visualization.quantity_support():
        fig, ax = plt.subplots(
            figsize=(7, 6),
            constrained_layout=True,
        )
        ax2 = ax.twiny()
        ax2.invert_xaxis()
        na.plt.plot(wavelength, depth, ax=ax, color="tab:green");
        na.plt.plot(energy, depth, ax=ax2, color="tab:green");
        ax.axhline(ccd.thickness_implant, linestyle="--", color="black")
        ax.text(
            x=0.1,
            y=ccd.thickness_implant / 10,
            s="PCC region",
            transform=ax.get_yaxis_transform(),
        )
        ax.set_xscale("log");
        ax2.set_xscale("log");
        ax.set_yscale("log");
        ax2.set_yscale("log");
        ax.set_xlabel(f"wavelength ({ax.get_xlabel()})");
        ax2.set_xlabel(f"energy ({ax2.get_xlabel()})", labelpad=16)
        ax.set_ylabel(f"penetration depth ({depth.unit:latex_inline})");

    path = default_path / "penetration_depth.svg"

    fig.savefig(path)

    return path
