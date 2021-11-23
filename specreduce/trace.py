# Licensed under a 3-clause BSD style license - see ../licenses/LICENSE.rst

from dataclasses import dataclass

# import numpy as np

# from scipy.signal import find_peaks

# from specreduce.core import SpecreduceOperation


@dataclass
class Aperture:
    """
    Class describing an aperture.

    Attributes
    ----------

    upper : float (default: 5.)
        Location of the upper edge of the aperture relative to the center

    lower : float (default: 5.)
        Location of the lower edge of the aperture relative to the center

    Properties
    ----------

    width : float
        Total width of the aperture
    """
    upper: float = 5.
    lower: float = -5.

    @property
    def width(self):
        return abs(self.upper - self.lower)

    @width.setter
    def width(self, newwidth):
        """
        Adjust aperture to a new width. Assumes aperture will be symmetric about the center.
        """
        if newwidth > 0:
            self.upper = 0.5 * newwidth
            self.lower = -self.upper
        else:
            raise ValueError(f"New width must be positive: {newwidth}")
