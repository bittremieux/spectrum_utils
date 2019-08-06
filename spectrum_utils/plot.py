import itertools
import math
from typing import Dict, Optional, Union

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw

from spectrum_utils.spectrum import MsmsSpectrum, MoleculeFragmentAnnotation, PeptideFragmentAnnotation


colors = {'a': '#388E3C', 'b': '#1976D2', 'c': '#00796B',
          'x': '#7B1FA2', 'y': '#D32F2F', 'z': '#F57C00',
          'unknown': '#212121', 'mol': '#212121', None: '#212121'}
zorders = {'a': 3, 'b': 4, 'c': 3, 'x': 3, 'y': 4, 'z': 3, 'unknown': 2,
           'mol': 5, None: 1}


_mol_annotation_size = 200


def _annotate_ion(mz, intensity, annotation, color_ions, annotate_ions,
                  annotation_kws, ax):
    # No annotation -> just return peak styling information.
    if annotation is None:
        return colors.get(None), zorders.get(None)
    # Else: add the textual or figure annotation.
    else:
        color = (colors.get(annotation.ion_type) if color_ions else
                 colors.get(None))
        zorder = zorders.get(annotation.ion_type)

        if annotate_ions:
            annotation_pos = intensity
            if annotation_pos > 0:
                annotation_pos += 0.02
            if type(annotation) == PeptideFragmentAnnotation:
                ax.text(mz, annotation_pos, str(annotation), color=color,
                        **annotation_kws)
            elif type(annotation) == MoleculeFragmentAnnotation:
                im = _smiles_to_im(annotation.smiles)
                x_range = ax.get_xlim()[1] - ax.get_xlim()[0]
                width = (x_range / 3) / _mol_annotation_size * im.shape[1]
                height = (1/3) / _mol_annotation_size * im.shape[0]
                ax.imshow(im, aspect='auto', extent=(mz - width // 2,
                                                     mz + width // 2,
                                                     annotation_pos,
                                                     annotation_pos + height),
                          zorder=zorder)

        return color, zorder


def _smiles_to_im(smiles):
    # Draw the molecule and make the white background transparent.
    # https://stackoverflow.com/a/54148416
    options = Draw.DrawingOptions()
    options.dotsPerAngstrom = 100
    im = Draw.MolToImage(Chem.MolFromSmiles(smiles),
                         size=(_mol_annotation_size, _mol_annotation_size),
                         options=options)
    im = np.asarray(im).copy()
    im[:, :, 3] = (255 * (im[:, :, :3] != 255).any(axis=2)).astype(np.uint8)
    # Crop the image by removing white lines.
    bottom, top, left, right = 0, im.shape[0] - 1, 0, im.shape[1] - 1
    while bottom < im.shape[0] and im[bottom:bottom+1, :, :3].min() == 255:
        bottom += 1
    while top > 0 and im[top-1:top, :, :3].min() == 255:
        top -= 1
    while left < im.shape[1] and im[:, left:left+1, :3].min() == 255:
        left += 1
    while right > 0 and im[:, right-1:right, :3].min() == 255:
        right -= 1
    return im[bottom:top:, left:right, :]


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
            mz, peak_intensity if not mirror_intensity else -peak_intensity,
            annotation, color_ions, annotate_ions, annotation_kws, ax)
        ax.plot([mz, mz], [0, peak_intensity], color=color, zorder=zorder)

    ax.xaxis.set_minor_locator(mticker.AutoLocator())
    ax.yaxis.set_minor_locator(mticker.AutoLocator())
    ax.xaxis.set_minor_locator(mticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator())
    if grid in (True, 'both', 'major'):
        ax.grid(b=True, which='major', color='#9E9E9E',
                linestyle='--', linewidth=1.0)
    if grid in (True, 'both', 'minor'):
        ax.grid(b=True, which='minor', color='#9E9E9E',
                linestyle='--', linewidth=0.5)
    ax.set_axisbelow(True)

    ax.tick_params(axis='both', which='both', labelsize='small')

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
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(
        lambda x, pos: f'{abs(x):.0%}'))

    return ax
