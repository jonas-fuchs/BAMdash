from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("bamdash")
except PackageNotFoundError:
    __version__ = "unknown"
