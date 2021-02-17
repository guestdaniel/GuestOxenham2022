"""
This script estimates tuning curves for simulated auditory nerve fibers at a range of CFs and saves the results
to disk for further analysis.
"""
import apcmodels.synthesis as sy
import apcmodels.simulation as si
import apcmodels.anf as anf
import numpy as np
from scipy.interpolate import interp1d
import config as cfg
import os


# Define function to estimate tuning curves
def estimate_tuning_curve(cf, levels, freqs, fs, model):
    """
    Calculates the mean firing rate of an auditory nerve model simulation responding to a short pure tone at a range of
    levels and frequencies centered around a cf.

    Parameters:
        cf (float): characteristic frequency to probe in Hz
        levels (list, ndarray): list or array of levels to test in dB SPL
        freqs (ndarray): list or array of frequencies to test, in octaves re: cf
        fs (int): sampling rate of the stimulus and model simulation in Hz
        model: model object from apcmodels.anf

    Returns:
        output (ndarray): array of mean firing rates of shape (n_freq, n_level)
    """
    # Write simple wrapper around runfunc to average all simulations over time as they are run
    def avg_wrapper(runfunc):
        def inner(params):
            return np.mean(runfunc(params))  # average output of runfunc over time
        return inner

    # Encode fixed stimulus params in dict
    params = si.Parameters(freq=cf, dur=0.10, dur_ramp=0.020, phase=0, fs=fs)  # encode fixed params
    params.wiggle('level', levels)  # wiggle levels
    params.wiggle('freq', cf*2**freqs)  # wiggle freqs

    # Add model params to dict
    params.append(['cf_low', 'cf_high', 'n_cf', 'fs'], [cf, cf, 1, fs])

    # Synthesize stimuli
    synth = sy.PureTone()
    stimuli = synth.synthesize_sequence(params)
    params.add_inputs(stimuli)

    # Run model
    sim = model()
    results = sim.run(params, runfunc=avg_wrapper(sim.simulate))

    # Return average of each response
    return results.astype('float32')  # convert to 32-bit precision floats


# Define parameters
fs = int(200e3)  # sampling rate in Hz
cfs = 10**np.linspace(np.log10(200), np.log10(20000), 25)  # CFs for which we will measure tuning curves
levels = np.linspace(-10, 40, num=25)  # range of levels over which we will estimate tuning curves
freqs = np.linspace(-1, 1, 25)  # range of frequencies over which we will estimate tuning curves\

# Save cfs, levels, and freqs to disk
np.save(os.path.join(cfg.root_directory, 'nofigure/tuning_curves/cfs.npy'), cfs)
np.save(os.path.join(cfg.root_directory, 'nofigure/tuning_curves/freqs.npy'), freqs)
np.save(os.path.join(cfg.root_directory, 'nofigure/tuning_curves/levels.npy'), levels)

# Loop over models and model names
for model, model_name in zip([anf.AuditoryNerveHeinz2001Numba], ['Heinz2001']):
    # Compute frequency-level profile for each cf
    tuning_curves = []
    for cf in cfs:
        tuning_curve = estimate_tuning_curve(cf, levels, freqs, fs, anf.AuditoryNerveHeinz2001Numba)
        tuning_curves.append(tuning_curve)
    # Save neural responses to disk
    np.save(os.path.join(cfg.root_directory, 'nofigure/tuning_curves/', model_name, '.npy'), tuning_curves)



