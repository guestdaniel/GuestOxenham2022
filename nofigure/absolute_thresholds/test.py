"""
This script estimates rate-level functions for simulated auditory nerve fibers at a range of CFs and saves the results
to disk for further analysis.
"""
import apcmodels.synthesis as sy
import apcmodels.simulation as si
import apcmodels.anf as anf
import numpy as np
from scipy.interpolate import interp1d
import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.getcwd())
import util
import os

print(os.getcwd())
