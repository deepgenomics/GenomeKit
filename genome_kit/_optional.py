from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version


def require_polars():
    """Import Polars if available, otherwise provide helpful error messages.

    Also checks for compability on MacOS with Apple Silicon, which may require
    an additional package if running Python under Rosetta translation.
    """
    try:
        import polars as pl

        if check_under_rosetta():
            if not check_rtcompat():
                raise ImportError(
                    "Polars is not compatible with Apple Silicon.\n"
                    "Please install with `pip install genomekit[df-mac]` to include "
                    "the polars-runtime-compat package required for Rosetta "
                    "translation."
                )
    except ModuleNotFoundError as e:
        raise ImportError(
            "Optional dependency 'polars' is required for this functionality. Please "
            "install with `pip install genomekit[df]`.\n"
            "If you are running this on MacOS with Apple Silicon, please install with "
            "`pip install genomekit[df-mac]` to include the polars-runtime-compat "
            "package required for Rosetta translation."
        ) from e

    return pl


def check_under_rosetta():
    """Check if program is running under Rosetta translation on Apple Silicon.

    The default version of Polars is incompatible with Rosetta, and requires
    polars-runtime-compat to be installed.

    Can be checked with the sysctl.proc_translated flag in sysctl.
    See https://developer.apple.com/documentation/apple-silicon/about-the-rosetta-translation-environment#Determine-Whether-Your-App-Is-Running-as-a-Translated-Binary
    """
    import subprocess

    try:
        result = subprocess.run(
            ["sysctl", "-n", "sysctl.proc_translated"],
            capture_output=True,
            text=True,
            check=True,
        )
        # output will be 0 if running commnad directly from the terminal, and 1 if
        # running through Python under Rosetta translation
        return result.stdout.strip() == "1"
    except subprocess.CalledProcessError:
        # sysctl.proc_translated won't exist on non-Apple Silicon machines
        return False


def check_rtcompat():
    """Check if polars-runtime-compat is installed.

    Required for Polars to run on MacOS machines under Rosetta translation.
    """
    try:
        version("polars-runtime-compat")
        return True
    except PackageNotFoundError:
        return False
