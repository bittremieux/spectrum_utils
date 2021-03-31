import itertools
import math
from typing import Dict, Optional, Tuple, Union

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

from spectrum_utils.spectrum import MsmsSpectrum, MoleculeFragmentAnnotation,\
    PeptideFragmentAnnotation


colors = {'a': '#388E3C', 'b': '#1976D2', 'c': '#00796B',
          'x': '#7B1FA2', 'y': '#D32F2F', 'z': '#F57C00',
          'unknown': '#212121', 'molecule': '#212121', None: '#212121'}
zorders = {'a': 3, 'b': 4, 'c': 3, 'x': 3, 'y': 4, 'z': 3, 'unknown': 2,
           'molecule': 5, None: 1}


def _annotate_ion(mz: float, intensity: float,
                  annotation: Optional[Union[MoleculeFragmentAnnotation,
                                             PeptideFragmentAnnotation]],
                  color_ions: bool, annotate_ions: bool,
                  annotation_kws: Dict[str, object], ax: plt.Axes)\
        -> Tuple[str, int]:
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

    Returns
    -------
    Tuple[str, int]
        A tuple of the annotation's color as a hex string and the annotation's
        zorder.
    """
    # No annotation -> Just return peak styling information.
    if annotation is None:
        return colors.get(None), zorders.get(None)
    # Else: Add the textual annotation.
    else:
        color = (colors.get(annotation.ion_type) if color_ions else
                 colors.get(None))
        zorder = zorders.get(annotation.ion_type)

        if annotate_ions:
            annotation_pos = intensity
            if annotation_pos > 0:
                annotation_pos += 0.02
            kws = annotation_kws.copy()
            del kws['zorder']
            ax.text(mz, annotation_pos, str(annotation), color=color,
                    zorder=zorder, **kws)

        return color, zorder


def spectrum(spec: MsmsSpectrum, color_ions: bool = True,
             annotate_ions: bool = True, annot_kws: Optional[Dict] = None,
             mirror_intensity: bool = False, grid: Union[bool, str] = True,
             ax: Optional[plt.Axes] = None) -> plt.Axes:
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

    min_mz = max(0, math.floor(spec.mz[0] / 100 - 1) * 100)
    max_mz = math.ceil(spec.mz[-1] / 100 + 1) * 100
    ax.set_xlim(min_mz, max_mz)
    ax.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1.))
    y_max = 1.15 if annotate_ions else 1.05
    ax.set_ylim(*(0, y_max) if not mirror_intensity else (-y_max, 0))

    max_intensity = spec.intensity.max()
    annotations = (spec.annotation if spec.annotation is not None else
                   itertools.repeat(None))
    annotation_kws = {
        'horizontalalignment': 'left' if not mirror_intensity else 'right',
        'verticalalignment': 'center', 'rotation': 90,
        'rotation_mode': 'anchor', 'zorder': 5}
    if annot_kws is not None:
        annotation_kws.update(annot_kws)
    for mz, intensity, annotation in zip(spec.mz, spec.intensity, annotations):
        peak_intensity = intensity / max_intensity
        if mirror_intensity:
            peak_intensity *= -1

        color, zorder = _annotate_ion(
            mz, peak_intensity, annotation, color_ions, annotate_ions,
            annotation_kws, ax)
        ax.plot([mz, mz], [0, peak_intensity], color=color, zorder=zorder)

    ax.xaxis.set_minor_locator(mticker.AutoLocator())
    ax.yaxis.set_minor_locator(mticker.AutoLocator())
    ax.xaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())
    if grid in (True, 'both', 'major'):
        ax.grid(b=True, which='major', color='#9E9E9E', linewidth=0.2)
    if grid in (True, 'both', 'minor'):
        ax.grid(b=True, which='minor', color='#9E9E9E', linewidth=0.2)
    ax.set_axisbelow(True)

    ax.tick_params(axis='both', which='both', labelsize='small')
    y_ticks = ax.get_yticks()
    ax.set_yticks(y_ticks[y_ticks <= 1.])

    ax.set_xlabel('m/z', style='italic')
    ax.set_ylabel('Intensity')

    return ax


def mirror(spec_top: MsmsSpectrum, spec_bottom: MsmsSpectrum,
           spectrum_kws: Optional[Dict] = None, ax: Optional[plt.Axes] = None)\
        -> plt.Axes:
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

    ax.axhline(0, color='#9E9E9E', zorder=10)

    # Update axes so that both spectra fit.
    min_mz = max([0, math.floor(spec_top.mz[0] / 100 - 1) * 100,
                  math.floor(spec_bottom.mz[0] / 100 - 1) * 100])
    max_mz = max([math.ceil(spec_top.mz[-1] / 100 + 1) * 100,
                  math.ceil(spec_bottom.mz[-1] / 100 + 1) * 100])
    ax.set_xlim(min_mz, max_mz)
    ax.yaxis.set_major_locator(mticker.AutoLocator())
    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(
        lambda x, pos: f'{abs(x):.0%}'))

    return ax
