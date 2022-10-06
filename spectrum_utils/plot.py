import functools
import itertools
import math
from typing import Callable, Dict, Iterable, Optional, Tuple, Union

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

import spectrum_utils.fragment_annotation as fa
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
    annotation: Optional[fa.FragmentAnnotation],
    color_ions: bool,
    annot_fmt: Optional[Callable],
    annot_kws: Dict[str, object],
    ax: plt.Axes,
) -> Tuple[str, int]:
    """
    Annotate a specific fragment peak.

    Parameters
    ----------
    mz : float
        The peak's m/z value (position of the annotation on the x axis).
    intensity : float
        The peak's intensity (position of the annotation on the y axis).
    annotation : Optional[fa.FragmentAnnotation]
        The annotation that will be plotted.
    color_ions : bool
        Flag whether to color the peak annotation or not.
    annot_fmt : Optional[Callable]
        Function to format the peak annotations. See `FragmentAnnotation` for
        supported elements. By default, only canonical b and y peptide fragments
        are annotated. If `None`, no peaks are annotated.
    annot_kws : Dict[str, object]
        Keyword arguments for `ax.text` to customize peak annotations.
    ax : plt.Axes
        Axes instance on which to plot the annotation.

    Returns
    -------
    Tuple[str, int]
        A tuple of the annotation's color as a hex string and the annotation's
        zorder.
    """
    ion_type = annotation.ion_type[0] if annotation is not None else None
    color = colors.get(ion_type if color_ions else None)
    zorder = zorders.get(ion_type)
    if annot_fmt is not None:
        y = intensity + 0.02 * (intensity > 0)
        kws = annot_kws.copy()
        kws.update(dict(color=color, zorder=zorder))
        ax.text(mz, y, annot_fmt(annotation), **kws)
    return color, zorder


def annotate_ion_type(
    annotation: fa.FragmentAnnotation, ion_types: Iterable[str]
) -> str:
    """
    Convert a `FragmentAnnotation` to a string for annotating peaks in a
    spectrum plot.

    This function will only annotate singly-charged, mono-isotopic canonical
    peaks with the given ion type(s).

    Parameters
    ----------
    annotation : fa.FragmentAnnotation
        The peak's fragment annotation.
    ion_types : Iterable[str]
        Accepted ion types to annotate.

    Returns
    -------
    str
        The peak's annotation string.
    """
    if (
        annotation.ion_type[0] in ion_types
        and annotation.neutral_loss is None
        and annotation.isotope == 0
        and annotation.charge == 1
    ):
        return annotation.ion_type
    else:
        return ""


def spectrum(
    spec: MsmsSpectrum,
    *,
    color_ions: bool = True,
    annot_fmt: Optional[Callable] = functools.partial(
        annotate_ion_type, ion_types="by"
    ),
    annot_kws: Optional[Dict] = None,
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
    annot_fmt : Optional[Callable]
        Function to format the peak annotations. See `FragmentAnnotation` for
        supported elements. By default, only canonical b and y peptide fragments
        are annotated. If `None`, no peaks are annotated.
    annot_kws : Optional[Dict], optional
        Keyword arguments for `ax.text` to customize peak annotations.
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

    ax.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1.0))
    ax.set_ylim(*(0, 1) if not mirror_intensity else (-1, 0))

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

    if len(spec.mz) == 0:
        return ax

    round_mz = 50
    max_mz = math.ceil(spec.mz[-1] / round_mz + 1) * round_mz
    ax.set_xlim(0, max_mz)

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
            annot_fmt,
            annotation_kws,
            ax,
        )
        ax.plot([mz, mz], [0, peak_intensity], color=color, zorder=zorder)

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

    max_mz_top = spec_top.mz[-1] if len(spec_top.mz) > 0 else 1
    max_mz_bottom = spec_bottom.mz[-1] if len(spec_bottom.mz) > 0 else 1
    # Update axes so that both spectra fit.
    round_mz = 50
    max_mz = max(
        [
            math.ceil(max_mz_top / round_mz + 1) * round_mz,
            math.ceil(max_mz_bottom / round_mz + 1) * round_mz,
        ]
    )
    ax.set_xlim(0, max_mz)
    ax.yaxis.set_major_locator(mticker.AutoLocator())
    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax.yaxis.set_major_formatter(
        mticker.FuncFormatter(lambda x, pos: f"{abs(x):.0%}")
    )

    return ax
