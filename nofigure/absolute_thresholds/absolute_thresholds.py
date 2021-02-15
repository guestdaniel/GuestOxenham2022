"""
This script estimates rate-level functions for simulated auditory nerve fibers at a range of CFs and saves the results
to disk for further analysis.
"""
import apcmodels.synthesis as sy
import apcmodels.simulation as si
import apcmodels.anf as anf
import numpy as np
from scipy.interpolate import interp1d
import config as cfg
import os


# Define function to estimate rate-level functions
def estimate_rate_level_function(cf, levels, fs, model):
    """
    Calculates the mean firing rate of an auditory nerve model simulation responding to a short pure tone at a range
    of levels at a given CF

    Parameters:
        cf (float): characteristic frequency to probe in Hz
        levels (list, ndarray): list or array of levels to test in dB SPL
        fs (int): sampling rate of the stimulus and model simulation in Hz
        model: model object from apcmodels.anf

    Returns:
        results (list): list of mean firing rates in spikes per second at each level in levels
    """
    # Encode fixed stimulus params in dict
    params = si.Parameters(freq=cf, dur=0.10, dur_ramp=0.020, phase=0, fs=fs)  # encode fixed params
    params.wiggle('level', levels)  # wiggle levels

    # Add model params to dict
    params.append(['cf_low', 'cf_high', 'n_cf', 'fs'], [cf, cf, 1, fs])

    # Synthesize stimuli
    synth = sy.PureTone()
    stimuli = synth.synthesize_sequence(params)
    params.add_inputs(stimuli)

    # Run model
    sim = model()
    results = sim.run(params)
    return [np.mean(result) for result in results]  # return mean of response at each level


# Parameters
fs = int(200e3)  # sampling rate in Hz
cfs = 10**np.linspace(np.log10(200), np.log10(20000), 25)  # CFs for which we will measure rate-level functions
levels = np.linspace(-10, 40, num=25)  # range of levels over which we will estimate rate-level functions

# Compute rate-level function for each CF and then estimate the level at which 1.05 x spontaneous rate is achieved
absolute_thresholds = []
for cf in cfs:
    rate_level_function = estimate_rate_level_function(cf, levels, fs, anf.AuditoryNerveHeinz2001Numba)
    absolute_thresholds.append(interp1d(rate_level_function, levels)(np.min([rate_level_function[0] * 1.05,
                                                                            np.max(rate_level_function)])))
# Save absolute thresholds to disk
np.save(os.path.join(cfg.root_directory, 'nofigure/absolute_thresholds/cfs.npy'), cfs)
np.save(os.path.join(cfg.root_directory, 'nofigure/absolute_thresholds/Heinz2001.npy'), absolute_thresholds)
