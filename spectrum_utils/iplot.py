import altair
import pandas as pd

from spectrum_utils.plot import colors
from spectrum_utils.spectrum import MsmsSpectrum


def spectrum(spec: MsmsSpectrum, color_ions: bool = True,
             annotate_ions: bool = True, mirror_intensity: bool = False,
             grid: bool = True) -> altair.LayerChart:
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
        annotation = (spec.annotation if annotate_ions else
                      [None] * len(spec.mz))
        color = [colors[a.ion_type if a is not None and color_ions else None]
                 for a in spec.annotation]
    else:
        annotation = [None] * len(spec.mz)
        color = [colors[None]] * len(spec.mz)
    calc_mz = [a.calc_mz if a is not None else None for a in annotation]
    fragment = [str(a) if a is not None else None for a in annotation]
    spec_df = pd.DataFrame({'exp_mz': spec.mz, 'intensity': intensity,
                            'calc_mz': calc_mz, 'fragment': fragment,
                            'color': color})

    x_axis = altair.X('exp_mz', axis=altair.Axis(title='m/z',
                                                 titleFontStyle='italic',
                                                 grid=grid),
                      scale=altair.Scale(nice=True, padding=5, zero=False))
    y_axis = altair.Y('intensity',
                      axis=altair.Axis(title='Intensity', format='%',
                                       grid=grid),
                      scale=altair.Scale(nice=True, padding=5))
    color = altair.Color('color', scale=None, legend=None)
    tooltip_not_annotated = [altair.Tooltip('exp_mz', format='.4f',
                                            title='Experimental m/z'),
                             altair.Tooltip('intensity', format='.0%',
                                            title='Intensity')]
    tooltip_annotated = [altair.Tooltip('fragment', title='Fragment'),
                         altair.Tooltip('exp_mz', format='.4f',
                                        title='Experimental m/z'),
                         altair.Tooltip('calc_mz', format='.4f',
                                        title='Theoretical m/z'),
                         altair.Tooltip('intensity', format='.0%',
                                        title='Intensity')]
    # Unannotated peaks.
    spec_plot = (altair.Chart(spec_df[spec_df['fragment'].isna()])
                 .mark_rule(size=2).encode(x=x_axis, y=y_axis, color=color,
                                           tooltip=tooltip_not_annotated))
    # Annotated peaks.
    spec_plot += (altair.Chart(spec_df[~spec_df['fragment'].isna()])
                  .mark_rule(size=2).encode(x=x_axis, y=y_axis, color=color,
                                            tooltip=tooltip_annotated))
    spec_plot += (altair.Chart(spec_df[~spec_df['fragment'].isna()])
                  .mark_text(align='right' if mirror_intensity else 'left',
                             angle=270, baseline='middle',
                             dx=-5 if mirror_intensity else 5)
                  .encode(x=x_axis, y=y_axis, text='fragment', color=color,
                          tooltip=tooltip_annotated))

    return spec_plot


def mirror(spec_top: MsmsSpectrum, spec_bottom: MsmsSpectrum,
           color_ions: bool = True, annotate_ions: bool = True)\
        -> altair.LayerChart:
    """
    Mirror plot two MS/MS spectra.

    Parameters
    ----------
    spec_top : MsmsSpectrum
        The spectrum to be plotted on the top.
    spec_bottom : MsmsSpectrum
        The spectrum to be plotted on the bottom.
    color_ions : bool, optional
        Flag indicating whether or not to color annotated fragment ions. The
        default is True.
    annotate_ions : bool, optional
        Flag indicating whether or not to annotate fragment ions. The default
        is True.

    Returns
    -------
    altair.LayerChart
        The Altair chart instance with the plotted spectrum.
    """
    # Top spectrum.
    spec_plot = spectrum(spec_top, color_ions, annotate_ions, False)
    # Mirrored bottom spectrum.
    spec_plot += spectrum(spec_bottom, color_ions, annotate_ions, True)

    spec_plot += (altair.Chart(pd.DataFrame({'sep': [0]}))
                  .mark_rule(size=3).encode(
                      y='sep', color=altair.value('lightGray')))

    return spec_plot
