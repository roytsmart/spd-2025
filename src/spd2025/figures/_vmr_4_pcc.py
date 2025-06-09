import pathlib
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.visualization
import named_arrays as na
import optika
import ccd_snr
from ._path import default_path


def vmr_pcc() -> pathlib.Path:

    ccd = ccd_snr.ccd()

    wavelength = na.geomspace(1.5, 10000, axis="wavelength", num=1001) << u.AA

    energy = wavelength.to(u.eV, equivalencies=u.spectral())

    rays = optika.rays.RayVectorArray(
        wavelength=wavelength,
        direction=na.Cartesian3dVectorArray(0, 0, 1),
    )

    normal = na.Cartesian3dVectorArray(0, 0, -1)

    iqy = ccd.quantum_yield_ideal(wavelength)
    cce = ccd.charge_collection_efficiency(rays, normal)
    eqe = ccd.quantum_efficiency_effective(rays, normal)
    qe = ccd.quantum_efficiency(rays, normal)

    absorbance = ccd.absorbance(rays, normal).average

    f = ccd.fano_factor(wavelength)

    n0 = ccd.cce_backsurface
    a = optika.chemicals.Chemical("Si").absorption(wavelength)
    W = ccd.thickness_implant
    aW = (a * W).to(u.dimensionless_unscaled).value

    mean_n = iqy
    var_n = f * mean_n
    mean_p = cce
    var_p = 2 * np.exp(-aW) * np.square((n0 - 1) / aW) * (np.sinh(aW) - aW)
    mean_p2 = np.square(mean_p)
    mean_n2 = np.square(mean_n)
    mean_i = mean_n * mean_p
    var_exp = (var_n * var_p) + (var_n * mean_p2) + (var_p * mean_n2)
    exp_var = mean_n * (mean_p - (var_p + mean_p2)) * u.electron / u.photon
    var_i = var_exp + exp_var

    vmr_simple = 1 / eqe * u.photon
    vmr_shot = 1 / absorbance * u.photon
    vmr_fano = f * u.photon / qe
    vmr_pcc = (var_i / mean_i * u.photon) / qe - vmr_fano

    color_simple = "gray"

    with astropy.visualization.quantity_support():
        fig, ax = plt.subplots(
            figsize=(7, 6),
            constrained_layout=True,
        )
        ax2 = ax.twiny()
        ax2.invert_xaxis()
        na.plt.plot(wavelength, vmr_simple, ax=ax, color=color_simple, label="simple")
        na.plt.plot(energy, vmr_simple, ax=ax2, color=color_simple)
        na.plt.plot(wavelength, vmr_shot, ax=ax, color="tab:blue", label="shot only")
        na.plt.plot(wavelength, vmr_fano, ax=ax, color="tab:green", label="Fano only")
        na.plt.plot(wavelength, vmr_pcc, ax=ax, color="tab:orange", label="PCC only")
        ax.set_xscale("log")
        ax2.set_xscale("log")
        ax.set_xlabel(f"wavelength ({ax.get_xlabel()})")
        ax2.set_xlabel(f"energy ({ax2.get_xlabel()})", labelpad=16)
        ax.set_ylabel(f"variance-to-mean ratio ({ax.get_ylabel()})")
        ax.set_ylim(bottom=-1, top=15)
        ax.legend()

    path = default_path / "vmr_pcc.svg"

    fig.savefig(path)

    return path
