import itertools
import math

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

from spectrum_utils.spectrum import MsmsSpectrum


colors = {'a': '#388E3C', 'b': '#1976D2', 'y': '#D32F2F',
          'unknown': '#212121', None: '#212121'}
zorders = {'a': 3, 'b': 4, 'c': 3, 'x': 3, 'y': 4, 'z': 3, 'unknown': 2,
           None: 1}


def spectrum(spec: MsmsSpectrum, color_ions: bool = True,
             annotate_ions: bool = True, ax: plt.Axes = None) -> plt.Axes:
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
    ax : plt.Axes, optional
        Axes instance on which to plot the spectrum. If None the current Axes
        instance is used.

    Returns
    -------
    plt.Axes
        The matplotlib Axes instance on which the spectrum is plotted.
    """
    if ax is None:
        ax = plt.gca()

    max_intensity = spec.intensity.max()
    annotations = (spec.annotation if spec.annotation is not None else
                   itertools.repeat(None))
    for mz, intensity, annotation in zip(spec.mz, spec.intensity, annotations):
        ion_type = annotation.ion_type if annotation is not None else None
        color = colors.get(ion_type) if color_ions else colors.get(None)
        zorder = zorders.get(ion_type)

        ax.plot([mz, mz], [0, intensity / max_intensity], color=color,
                zorder=zorder)

        if annotate_ions and annotation is not None:
            ax.text(mz + 5, intensity / max_intensity + 0.02, str(annotation),
                    color=color, zorder=5, rotation=90, rotation_mode='anchor')

    min_mz = max(0, math.floor(spec.mz[0] / 100 - 1) * 100)
    max_mz = math.ceil(spec.mz[-1] / 100 + 1) * 100
    ax.set_xlim(min_mz, max_mz)
    ax.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1.))
    ax.set_ylim(0, 1.15 if annotate_ions else 1.05)

    ax.xaxis.set_minor_locator(mticker.AutoLocator())
    ax.yaxis.set_minor_locator(mticker.AutoLocator())
    ax.xaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax.grid(b=True, which='major', color='#9E9E9E',
            linestyle='--', linewidth=1.0)
    ax.grid(b=True, which='minor', color='#9E9E9E',
            linestyle='--', linewidth=0.5)
    ax.set_axisbelow(True)

    ax.tick_params(axis='both', which='both', labelsize='small')

    ax.set_xlabel('m/z')
    ax.set_ylabel('Intensity')

    return ax
