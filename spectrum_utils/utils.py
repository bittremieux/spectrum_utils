from typing import Union

import numba as nb
import numpy as np


@nb.njit(fastmath=True, cache=True)
def mass_diff(mz1, mz2, mode_is_da):
    """
    Calculate the mass difference(s).

    Parameters
    ----------
    mz1
        First m/z value(s).
    mz2
        Second m/z value(s).
    mode_is_da : bool
        Mass difference in Dalton (True) or in ppm (False).

    Returns
    -------
        The mass difference(s) between the given m/z values.
    """
    return mz1 - mz2 if mode_is_da else (mz1 - mz2) / mz2 * 10**6


def da_to_ppm(
    delta_mz: Union[int, np.ndarray], mz: Union[int, np.ndarray]
) -> Union[int, np.ndarray]:
    """
    Convert a mass difference in Dalton to ppm.

    Parameters
    ----------
    delta_mz : int or np.ndarray
        Mass difference in Dalton.
    mz : int or np.ndarray
        m/z value of peak.

    Returns
    -------
    int or np.ndarray

    """
    return delta_mz / mz * 1e6


def ppm_to_da(
    delta_mz: Union[int, np.ndarray], mz: Union[int, np.ndarray]
) -> Union[int, np.ndarray]:
    """
    Convert a mass difference in ppm to Dalton.

    Parameters
    ----------
    delta_mz : int or np.ndarray
        Mass difference in ppm.
    mz : int or np.ndarray
        m/z value of peak.

    Returns
    -------
    int or np.ndarray

    """
    return delta_mz / 1e6 * mz
