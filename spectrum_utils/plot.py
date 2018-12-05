import math

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker


colors = {'a': '#388E3C', 'b': '#1976D2', 'y': '#D32F2F',
          'unknown': '#212121', None: '#212121'}
zorders = {'a': 3, 'b': 4, 'c': 3, 'x': 3, 'y': 4, 'z': 3, 'unknown': 2,
           None: 1}


def spectrum(spec, color_ions=True, annotate_ions=True, ax=None):
    if ax is None:
        ax = plt.gca()

    max_intensity = spec.intensity.max()
    for mz, intensity, annotation in zip(spec.mz, spec.intensity,
                                         spec.annotation):
        ion_type = annotation.ion_type if annotation is not None else None
        color = colors.get(ion_type) if color_ions else colors.get(None)
        zorder = zorders.get(ion_type)

        ax.plot([mz, mz], [0, intensity / max_intensity], color=color,
                zorder=zorder)

        if annotate_ions and annotation is not None:
            ax.text(mz + 5, intensity / max_intensity + 0.02, str(annotation),
                    color=color, zorder=5, rotation=90, rotation_mode='anchor')

    if spec.peptide is not None:
        title = f'{spec.peptide}'
        if spec.precursor_charge is not None:
            title += f'/{spec.precursor_charge}'
        ax.text(0.5, 1.06, title,
                horizontalalignment='center', verticalalignment='bottom',
                fontsize='x-large', fontweight='bold', transform=ax.transAxes)
    ax.text(0.5, 1.02,
            f'Precursor m/z: {spec.precursor_mz:.4f}, '
            f'Charge: {spec.precursor_charge}',
            horizontalalignment='center', verticalalignment='bottom',
            fontsize='large', transform=ax.transAxes)

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
