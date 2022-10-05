import itertools
import math
from typing import Dict, Optional, Tuple, Union

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

from spectrum_utils.fragment_annotation import FragmentAnnotation
from spectrum_utils.spectrum import MsmsSpectrum


colors = {
    "a": "#388E3C",
    "b": "#1976D2",
    "c": "#00796B",
    "x": "#7B1FA2",
    "y": "#D32F2F",
    "z": "#F57C00",
    "m": "#FBC02D",
    "I": "#455A64",
    "p": "#512DA8",
    "?": "#212121",
    "f": "#212121",
    None: "#212121",
}
zorders = {
    "a": 3,
    "b": 4,
    "c": 3,
    "x": 3,
    "y": 4,
    "z": 3,
    "m": 2,
    "I": 3,
    "p": 3,
    "?": 2,
    "f": 5,
    None: 1,
}


def _annotate_ion(
    mz: float,
    intensity: float,
    annotation: Optional[FragmentAnnotation],
    color_ions: bool,
    annotate_ions: bool,
    annotation_kws: Dict[str, object],
    ax: plt.Axes,
    annot_fmt: Optional[str] = None,
) -> Tuple[str, int]:
    """
    Annotate a specific fragment peak.

    Parameters
    ----------
    mz : float
        The peak's m/z value (position of the annotation on the x axis).
    intensity : float
        The peak's intensity (position of the annotation on the y axis).
    annotation : Optional[MoleculeFragmentAnnotation,
                          PeptideFragmentAnnotation]
        The annotation that will be plotted.
    color_ions : bool
        Flag whether to color the peak annotation or not.
    annotate_ions : bool
        Flag whether to annotation the peak or not.
    annotation_kws : Dict
        Keyword arguments for `ax.text` to customize peak annotations.
    ax : plt.Axes
        Axes instance on which to plot the annotation.
    annot_fmt : Optional[str]
        Formatting string for peak annotations. Supported elements are
        '{ion_type}", '{neutral_loss}", '{isotope}", '{charge}", '{adduct}",
        '{analyte_number}", and '{mz_delta}". Example:
        "{ion_type}{neutral_loss}^{charge}".

    Returns
    -------
    Tuple[str, int]
        A tuple of the annotation's color as a hex string and the annotation's
        zorder.
    """
    ion_type = annotation.ion_type[0] if annotation is not None else None
    color = colors.get(ion_type if color_ions else None)
    zorder = zorders.get(ion_type)
    if annotate_ions and ion_type is not None and ion_type != "?":
        y = intensity + 0.02 * (intensity > 0)
        if annot_fmt is None:
            annot_str = str(annotation)
        else:
            annot_str = annot_fmt.format(
                ion_type=annotation.ion_type,
                neutral_loss=annotation.neutral_loss if annotation.neutral_loss is not None else "",
                isotope=annotation.isotope,
                charge=annotation.charge,
                adduct=annotation.adduct,
                analyte_number=annotation.analyte_number,
                mz_delta="".join(map(str, annotation.mz_delta)),
            )
        kws = annotation_kws.copy()
        kws.update(dict(color=color, zorder=zorder))
        ax.text(mz, y, annot_str, **kws)

    return color, zorder


def spectrum(
    spec: MsmsSpectrum,
    color_ions: bool = True,
    annotate_ions: bool = True,
    annot_kws: Optional[Dict] = None,
    annot_fmt: Optional[str] = None,
    mirror_intensity: bool = False,
    grid: Union[bool, str] = True,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """
    Plot an MS/MS spectrum.

    Parameters
    ----------
    spec : MsmsSpectrum
        The spectrum to be plotted.
    color_ions : bool, optional
        Flag indicating whether or not to color annotated fragment ions. The
        default is True.
    annotate_ions : bool, optional
        Flag indicating whether or not to annotate fragment ions. The default
        is True.
    annot_kws : Optional[Dict], optional
        Keyword arguments for `ax.text` to customize peak annotations.
    annot_fmt : Optional[str]
        Formatting string for peak annotations. Supported elements are
        '{ion_type}", '{neutral_loss}", '{isotope}", '{charge}", '{adduct}",
        '{analyte_number}", and '{mz_delta}". Example:
        "{ion_type}{neutral_loss}^{charge}".
    mirror_intensity : bool, optional
        Flag indicating whether to flip the intensity axis or not.
    grid : Union[bool, str], optional
        Draw grid lines or not. Either a boolean to enable/disable both major
        and minor grid lines or 'major'/'minor' to enable major or minor grid
        lines respectively.
    ax : Optional[plt.Axes], optional
        Axes instance on which to plot the spectrum. If None the current Axes
        instance is used.

    Returns
    -------
    plt.Axes
        The matplotlib Axes instance on which the spectrum is plotted.
    """
    if ax is None:
        ax = plt.gca()

    round_mz = 50
    max_mz = math.ceil(spec.mz[-1] / round_mz + 1) * round_mz
    ax.set_xlim(0, max_mz)
    ax.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1.0))
    ax.set_ylim(*(0, 1) if not mirror_intensity else (-1, 0))

    max_intensity = spec.intensity.max()
    annotations = (
        spec.annotation
        if spec.annotation is not None
        else itertools.repeat(None)
    )
    annotation_kws = {
        "horizontalalignment": "left" if not mirror_intensity else "right",
        "verticalalignment": "center",
        "rotation": 90,
        "rotation_mode": "anchor",
        "zorder": 5,
    }
    if annot_kws is not None:
        annotation_kws.update(annot_kws)
    for mz, intensity, annotation in zip(spec.mz, spec.intensity, annotations):
        peak_intensity = intensity / max_intensity
        if mirror_intensity:
            peak_intensity *= -1

        color, zorder = _annotate_ion(
            mz,
            peak_intensity,
            # Use the first annotation in case there are multiple options.
            annotation[0],
            color_ions,
            annotate_ions,
            annotation_kws,
            ax,
            annot_fmt,
        )
        ax.plot([mz, mz], [0, peak_intensity], color=color, zorder=zorder)

    ax.xaxis.set_minor_locator(mticker.AutoLocator())
    ax.yaxis.set_minor_locator(mticker.AutoLocator())
    ax.xaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())
    if grid in (True, "both", "major"):
        ax.grid(b=True, which="major", color="#9E9E9E", linewidth=0.2)
    if grid in (True, "both", "minor"):
        ax.grid(b=True, which="minor", color="#9E9E9E", linewidth=0.2)
    ax.set_axisbelow(True)

    ax.tick_params(axis="both", which="both", labelsize="small")

    ax.set_xlabel("m/z", style="italic")
    ax.set_ylabel("Intensity")

    return ax


def mirror(
    spec_top: MsmsSpectrum,
    spec_bottom: MsmsSpectrum,
    spectrum_kws: Optional[Dict] = None,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """
    Mirror plot two MS/MS spectra.

    Parameters
    ----------
    spec_top : MsmsSpectrum
        The spectrum to be plotted on the top.
    spec_bottom : MsmsSpectrum
        The spectrum to be plotted on the bottom.
    spectrum_kws : Optional[Dict], optional
        Keyword arguments for `plot.spectrum`.
    ax : Optional[plt.Axes], optional
        Axes instance on which to plot the spectrum. If None the current Axes
        instance is used.

    Returns
    -------
    plt.Axes
        The matplotlib Axes instance on which the spectra are plotted.
    """
    if ax is None:
        ax = plt.gca()

    if spectrum_kws is None:
        spectrum_kws = {}
    # Top spectrum.
    spectrum(spec_top, mirror_intensity=False, ax=ax, **spectrum_kws)
    y_max = ax.get_ylim()[1]
    # Mirrored bottom spectrum.
    spectrum(spec_bottom, mirror_intensity=True, ax=ax, **spectrum_kws)
    y_min = ax.get_ylim()[0]
    ax.set_ylim(y_min, y_max)

    ax.axhline(0, color="#9E9E9E", zorder=10)

    # Update axes so that both spectra fit.
    round_mz = 50
    max_mz = max(
        [
            math.ceil(spec_top.mz[-1] / round_mz + 1) * round_mz,
            math.ceil(spec_bottom.mz[-1] / round_mz + 1) * round_mz,
        ]
    )
    ax.set_xlim(0, max_mz)
    ax.yaxis.set_major_locator(mticker.AutoLocator())
    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax.yaxis.set_major_formatter(
        mticker.FuncFormatter(lambda x, pos: f"{abs(x):.0%}")
    )

    return ax
