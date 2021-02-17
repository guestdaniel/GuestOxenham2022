"""
This script takes tuning curves and extracts measures of q10 for each auditory nerve model.
"""
from scipy.interpolate import interp1d
import numpy as np
import config as cfg
import os


def estimate_q10(freqs, levels):
    """
    Takes a tuning curve and estimates sharpness in terms of q10

    Parameters:
        freqs (np.ndarray): an array of frequencies in Hz
        levels (np.ndarray): an array of levels in dB SPL

    Returns:
         q10 (float): sharpness in q10
    """
    f_peak = freqs[np.argmin(levels)]
    f_low = interp1d(levels[0:(np.argmin(levels) + 1)], freqs[0:(np.argmin(levels) + 1)])(np.min(levels) + 10)
    f_high = interp1d(levels[(np.argmin(levels) + 1):], freqs[(np.argmin(levels) + 1):])(np.min(levels) + 10)
    return f_peak / (f_high - f_low)


# Create storage and loop through and estimate q10
q10s = list()
for model_name in ['Heinz2001']:
    tuning_curves = np.load(os.path.join(cfg.root_directory, 'nofigure/tuning_curves/', model_name + '_tuning_curves' +
                                         '.npy'), allow_pickle=True)
    for tc in tuning_curves:
        q10s.append(estimate_q10(tc[0], tc[1]))
    # Save q10 estimates to disk
    np.save(os.path.join(cfg.root_directory, 'nofigure/tuning_curves/', model_name + '_q10s' + '.npy'),
            q10s)
