from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("bamdash")
except PackageNotFoundError:
    __version__ = "unknown"
