import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import named_arrays as na
import ccd_snr
from ._path import default_path

__all__ = [
    "differential_cce",
]


def differential_cce() -> list:

    thickness_oxide = 1 * u.um
    thickness_epitaxial = 16 * u.um
    thickness_pcc = 2.6 * u.um
    thickness_depletion = 6 * u.um

    cce_backsurface = ccd_snr.ccd().cce_backsurface

    unit = u.um
    thickness_oxide = thickness_oxide.to_value(unit)
    thickness_epitaxial = thickness_epitaxial.to_value(unit)
    thickness_pcc = thickness_pcc.to_value(unit)
    thickness_depletion = thickness_depletion.to_value(unit)
    thickness_fieldfree = thickness_epitaxial - thickness_depletion

    z = na.linspace(0, thickness_epitaxial, axis="z", num=1001)
    potential_recombination = np.exp(-((z - thickness_pcc) ** 2) / 2 / (1 / 2) ** 2)
    potential_depletion = np.minimum(
        1,
        -(1 / thickness_depletion) * (z - thickness_epitaxial),
    )
    potential = potential_recombination + potential_depletion

    cce = np.minimum(
        1,
        cce_backsurface + (1 - cce_backsurface) * z / thickness_pcc,
    )

    probability = np.exp(-z / thickness_pcc)

    fig, ax = plt.subplots(
        figsize=(8, 6),
        constrained_layout=True,
    )

    ax.margins(0.01)

    transform = ax.get_xaxis_transform()

    ymax = 0.7
    ax.axvspan(
        xmin=-thickness_oxide,
        xmax=0,
        ymax=ymax,
        facecolor="gray",
        edgecolor="black",
        label=r"SiO$_2$",
        zorder=5,
    )
    ax.axvspan(
        xmin=0,
        xmax=thickness_epitaxial,
        ymax=ymax,
        facecolor="lightgray",
        edgecolor="black",
        label=r"Si",
    )
    ax.axvspan(
        xmin=0,
        xmax=thickness_pcc,
        ymax=ymax,
        facecolor="none",
        edgecolor="gray",
        hatch="//",
    )
    # ax.axvspan(
    #     xmin=thickness_epitaxial - thickness_depletion,
    #     xmax=thickness_epitaxial,
    #     ymax=ymax,
    #     facecolor="none",
    #     edgecolor="gray",
    #     hatch="o",
    # )
    ax.axvspan(
        xmin=0,
        xmax=thickness_epitaxial,
        ymax=ymax,
        facecolor="none",
        edgecolor="black",
        zorder=10,
    )

    # color_potential = "tab:blue"
    # na.plt.plot(
    #     z,
    #     potential,
    #     ax=ax,
    #     zorder=10,
    #     color=color_potential,
    # )
    # ax.text(
    #     x=thickness_pcc + (thickness_fieldfree - thickness_pcc) / 2,
    #     y=1.03,
    #     s="electric\npotential",
    #     ha="center",
    #     va="bottom",
    #     color=color_potential,
    # )

    color_cce = "tab:orange"
    ax2 = ax.twinx()
    na.plt.plot(
        z,
        cce,
        ax=ax2,
        color=color_cce,
        zorder=10,
    )
    ax2.text(
        x=thickness_pcc + (thickness_fieldfree - thickness_pcc) / 2,
        y=0.98,
        s=r"$\eta(z)$",
        ha="center",
        va="top",
        color=color_cce,
    )

    color_probability = "tab:red"
    na.plt.plot(
        z,
        probability,
        ax=ax2,
        color=color_probability,
        zorder=5,
    )
    ax2.text(
        x=z[dict(z=300)].ndarray,
        y=probability[dict(z=300)].ndarray,
        s=" absorption probability",
        ha="left",
        va="bottom",
        color=color_probability,
    )

    kwargs_arrow = dict(
        dy=0,
        head_width=0.01,
        head_length=0.3,
        transform=transform,
        facecolor="black",
        edgecolor="none",
        # overhang=1,
        length_includes_head=True,
        linewidth=2,
    )

    bbox = dict(facecolor="white", edgecolor="none", boxstyle="square,pad=0.1")

    color_vline = "lightgray"

    x_D = thickness_epitaxial / 2
    y_D = 0.95
    ax.text(
        x=x_D,
        y=y_D,
        s="$D$",
        ha="center",
        va="center",
        transform=transform,
        bbox=bbox,
    )
    ax.arrow(
        x=x_D,
        y=y_D,
        dx=-x_D,
        **kwargs_arrow,
    )
    ax.arrow(
        x=x_D,
        y=y_D,
        dx=thickness_epitaxial - x_D,
        **kwargs_arrow,
    )
    ax.axvline(
        x=0,
        ymax=y_D + 0.02,
        color=color_vline,
        zorder=0,
    )
    ax.axvline(
        x=thickness_epitaxial,
        ymax=y_D + 0.02,
        color=color_vline,
        zorder=0,
    )

    x_delta = -thickness_oxide / 2
    y_delta = 0.9
    ax.text(
        x=x_delta,
        y=y_delta,
        s=r"$\delta$",
        ha="center",
        va="center",
        transform=transform,
        bbox=bbox,
        zorder=6,
    )
    ax.arrow(
        x=-thickness_oxide - 1,
        y=y_delta,
        dx=1,
        **kwargs_arrow,
    )
    ax.arrow(
        x=1,
        y=y_delta,
        dx=-1,
        **kwargs_arrow,
    )
    ax.axvline(
        x=-thickness_oxide,
        ymax=y_D + 0.02,
        color=color_vline,
        zorder=0,
    )

    x_W = thickness_pcc / 2
    y_W = 0.8
    ax.text(
        x=x_W,
        y=y_W,
        s="$W$",
        ha="center",
        va="center",
        transform=transform,
        bbox=bbox,
    )
    ax.arrow(
        x=x_W,
        y=y_W,
        dx=-thickness_pcc / 2,
        **kwargs_arrow,
    )
    ax.arrow(
        x=x_W,
        y=y_W,
        dx=thickness_pcc / 2,
        **kwargs_arrow,
    )
    ax.axvline(
        x=thickness_pcc,
        ymax=y_W + 0.02,
        color=color_vline,
        zorder=0,
    )

    # x_zd = thickness_epitaxial - thickness_depletion / 2
    # y_zd = 0.85
    # ax.text(
    #     x=x_zd,
    #     y=y_zd,
    #     s="$z_d$",
    #     ha="center",
    #     va="center",
    #     transform=transform,
    #     bbox=bbox,
    # )
    # ax.arrow(
    #     x=x_zd,
    #     y=y_zd,
    #     dx=-thickness_depletion / 2,
    #     **kwargs_arrow,
    # )
    # ax.arrow(
    #     x=x_zd,
    #     y=y_zd,
    #     dx=thickness_depletion / 2,
    #     **kwargs_arrow,
    # )
    # ax.axvline(
    #     x=thickness_epitaxial - thickness_depletion,
    #     ymax=y_zd + 0.02,
    #     color=color_vline,
    #     zorder=0,
    # )
    #
    # x_zf = thickness_fieldfree / 2
    # y_zf = 0.85
    # ax.text(
    #     x=x_zf,
    #     y=y_zf,
    #     s="$z_f$",
    #     ha="center",
    #     va="center",
    #     transform=transform,
    #     bbox=bbox,
    # )
    # ax.arrow(
    #     x=x_zf,
    #     y=y_zf,
    #     dx=-x_zf,
    #     **kwargs_arrow,
    # )
    # ax.arrow(
    #     x=x_zf,
    #     y=y_zf,
    #     dx=thickness_fieldfree - x_zf,
    #     **kwargs_arrow,
    # )

    ax.arrow(
        x=-6,
        y=0.35,
        dx=3,
        dy=0,
        head_width=0.04,
        head_length=1,
        transform=transform,
        facecolor="black",
        linewidth=0.5,
    )
    ax.text(
        x=-4,
        y=0.4,
        s="incident\nlight",
        va="bottom",
        ha="center",
        transform=transform,
    )

    ax.text(
        x=thickness_pcc / 2,
        y=0.03,
        s="PCC\nregion",
        transform=transform,
        ha="center",
        bbox=dict(facecolor="lightgray", edgecolor="none", boxstyle="square,pad=0.1"),
    )
    # ax.text(
    #     x=thickness_fieldfree + thickness_depletion / 2,
    #     y=0.5,
    #     s="depletion\nregion",
    #     transform=transform,
    #     ha="center",
    #     bbox=dict(facecolor="lightgray", edgecolor="none", boxstyle="square,pad=0.1"),
    # )
    # ax.text(
    #     x=thickness_pcc + (thickness_fieldfree - thickness_pcc) / 2,
    #     y=0.1,
    #     s="field-free\nregion",
    #     transform=transform,
    #     ha="center",
    #     bbox=dict(facecolor="lightgray", edgecolor="none", boxstyle="square,pad=0.1"),
    # )

    ax.legend(loc="lower left", handlelength=1)

    ax.set_xlabel(f"$z$ ({u.um:latex_inline})")
    # ax.set_ylabel("electric potential (arb. units)", color=color_potential)
    ax2.set_ylabel("differential CCE", color=color_cce)

    # ax.tick_params(axis="y", labelcolor=color_potential)
    ax2.tick_params(axis="y", labelcolor=color_cce)

    ax.set_ylim(top=3)
    ax2.set_ylim(0, 1.5)

    ax.set_axis_off()
    ax2.set_axis_off()

    path = default_path / "differential_cce.svg"

    fig.savefig(path)

    return path
