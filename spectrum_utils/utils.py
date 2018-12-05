def mass_diff(mz1, mz2, mode: str = 'Da'):
    """
    Calculate the mass difference(s).

    Parameters
    ----------
    mz1
        First m/z value(s).
    mz2
        Second m/z value(s).
    mode : {'Da', 'ppm'}
        Mass difference unit. Either 'Da' or 'ppm'.

    Returns
    -------
        The mass difference(s) between the given m/z values.
    """
    if mode == 'Da':
        return mz1 - mz2
    elif mode == 'ppm':
        return (mz1 - mz2) / mz2 * 10**6
