"""Font metrics for label layout.

Labels are typeset with Arial-compatible metrics (ACS 1996 document settings use Arial 10 pt).
Liberation Sans is metric-compatible with Arial, so either font gives identical advance widths.

Lookup order: the VISUALATLAS_FONT environment variable, then common Liberation Sans paths
(Linux), then the Arial that ships with macOS and Windows.
"""

import os
from functools import lru_cache

from fontTools.ttLib import TTFont

CANDIDATES = (
    "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",  # Debian / Ubuntu
    "/usr/share/fonts/liberation-sans/LiberationSans-Regular.ttf",  # Fedora
    "/usr/share/fonts/liberation/LiberationSans-Regular.ttf",  # Arch
    "/opt/homebrew/share/fonts/LiberationSans-Regular.ttf",
    "/System/Library/Fonts/Supplemental/Arial.ttf",  # macOS
    "/Library/Fonts/Arial.ttf",
    "C:/Windows/Fonts/arial.ttf",
)


def font_path():
    env = os.environ.get("VISUALATLAS_FONT")
    if env:
        if not os.path.exists(env):
            raise FileNotFoundError(f"VISUALATLAS_FONT points to a missing file: {env}")
        return env
    for p in CANDIDATES:
        if os.path.exists(p):
            return p
    raise FileNotFoundError(
        "No Arial-compatible font found. Install Liberation Sans (e.g. `apt install fonts-liberation`) "
        "or set VISUALATLAS_FONT to the path of a .ttf file."
    )


@lru_cache(maxsize=1)
def _metrics():
    f = TTFont(font_path())
    return f["head"].unitsPerEm, f.getBestCmap(), f["hmtx"]


def adv(c):
    """Advance width of character c, in units of the font size."""
    upm, cmap, hmtx = _metrics()
    g = cmap.get(ord(c))
    return hmtx[g][0] / upm if g else 0.6
