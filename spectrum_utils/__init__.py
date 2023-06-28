try:
    from importlib.metadata import version, PackageNotFoundError

    try:
        __version__ = version("spectrum_utils")
    except PackageNotFoundError:
        pass
except ImportError:
    from pkg_resources import get_distribution, DistributionNotFound

    try:
        __version__ = get_distribution("spectrum_utils").version
    except DistributionNotFound:
        pass


__all__ = [
    "fragment_annotation",
    "iplot",
    "plot",
    "proforma",
    "spectrum",
    "utils",
]
