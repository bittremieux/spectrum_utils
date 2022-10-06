import functools
import operator
from typing import Callable, Dict, Optional

try:
    import altair
    import pandas as pd
except ImportError:
    raise ImportError(
        "Missing dependencies for interactive plotting. Install using `pip "
        "install spectrum_utils[iplot]`, manually install Altair and Pandas, or"
        " use the default Matplotlib (`spectrum_utils.plot`) plotting backend."
    )

from spectrum_utils.plot import annotate_ion_type, colors
from spectrum_utils.spectrum import MsmsSpectrum


def spectrum(
    spec: MsmsSpectrum,
    *_,
    color_ions: bool = True,
    annot_fmt: Optional[Callable] = functools.partial(
        annotate_ion_type, ion_types="by"
    ),
    annot_kws: Optional[Dict] = None,
    mirror_intensity: bool = False,
    grid: bool = True,
) -> altair.LayerChart:
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
        Keyword arguments for `altair.Chart.mark_text` to customize peak
        annotations.
    mirror_intensity : bool, optional
        Flag indicating whether to flip the intensity axis or not.
    grid : bool, optional
        Draw grid lines or not.

    Returns
    -------
    altair.LayerChart
        The Altair chart instance with the plotted spectrum.
    """
    intensity = spec.intensity / spec.intensity.max()
    if mirror_intensity:
        intensity *= -1
    if spec.annotation is not None:
        annotations = list(map(operator.itemgetter(0), spec.annotation))
        peak_labels = map(annot_fmt, annotations)
        peak_colors = [
            colors.get(a.ion_type[0] if color_ions else None)
            for a in annotations
        ]
        mz_delta = [
            None if a.mz_delta is None else "".join(map(str, a.mz_delta))
            for a in annotations
        ]
        spec_df = pd.DataFrame(
            {
                "mz": spec.mz,
                "intensity": intensity,
                "fragment": peak_labels,
                "mz_delta": mz_delta,
                "color": peak_colors,
            }
        )
    else:
        spec_df = pd.DataFrame(
            {
                "mz": spec.mz,
                "intensity": intensity,
                "color": [colors[None]] * len(spec.mz),
            }
        )

    x_axis = altair.X(
        "mz",
        axis=altair.Axis(title="m/z", titleFontStyle="italic", grid=grid),
        scale=altair.Scale(nice=True, padding=5),
    )
    y_axis = altair.Y(
        "intensity",
        axis=altair.Axis(title="Intensity", format="%", grid=grid),
        scale=altair.Scale(nice=True),
    )
    color = altair.Color("color", scale=None, legend=None)
    tooltip_not_annotated = [
        altair.Tooltip("mz", format=".4f", title="m/z"),
        altair.Tooltip("intensity", format=".1%", title="Intensity"),
    ]
    tooltip_annotated = [
        altair.Tooltip("mz", format=".4f", title="m/z"),
        altair.Tooltip("intensity", format=".1%", title="Intensity"),
        altair.Tooltip("fragment", title="Fragment"),
        altair.Tooltip("mz_delta", title="m/z deviation"),
    ]
    # Unannotated peaks.
    mask_unannotated = spec_df["fragment"] == ""
    spec_plot = (
        altair.Chart(spec_df[mask_unannotated])
        .mark_rule(size=2)
        .encode(x=x_axis, y=y_axis, color=color, tooltip=tooltip_not_annotated)
    )
    # Annotated peaks.
    annotation_kws = {
        "align": "left" if not mirror_intensity else "right",
        "angle": 270,
        "baseline": "middle",
    }
    if annot_kws is not None:
        annotation_kws.update(annot_kws)
    spec_plot += (
        altair.Chart(spec_df[~mask_unannotated])
        .mark_rule(size=2)
        .encode(x=x_axis, y=y_axis, color=color, tooltip=tooltip_annotated)
    )
    spec_plot += (
        altair.Chart(spec_df[~mask_unannotated])
        .mark_text(dx=-5 if mirror_intensity else 5, **annotation_kws)
        .encode(
            x=x_axis,
            y=y_axis,
            text="fragment",
            color=color,
            tooltip=tooltip_annotated,
        )
    )

    return spec_plot


def mirror(
    spec_top: MsmsSpectrum,
    spec_bottom: MsmsSpectrum,
    spectrum_kws: Optional[Dict] = None,
    *_,
) -> altair.LayerChart:
    """
    Mirror plot two MS/MS spectra.

    Parameters
    ----------
    spec_top : MsmsSpectrum
        The spectrum to be plotted on the top.
    spec_bottom : MsmsSpectrum
        The spectrum to be plotted on the bottom.
    spectrum_kws : Optional[Dict], optional
        Keyword arguments for `iplot.spectrum`.
    *_
        Ignored, for consistency with the `plot.mirror` API.

    Returns
    -------
    altair.LayerChart
        The Altair chart instance with the plotted spectrum.
    """
    if spectrum_kws is None:
        spectrum_kws = {}
    # Top spectrum.
    spec_plot = spectrum(spec_top, mirror_intensity=False, **spectrum_kws)
    # Mirrored bottom spectrum.
    spec_plot += spectrum(spec_bottom, mirror_intensity=True, **spectrum_kws)

    spec_plot += (
        altair.Chart(pd.DataFrame({"sep": [0]}))
        .mark_rule(size=3)
        .encode(y="sep", color=altair.value("lightGray"))
    )

    return spec_plot
