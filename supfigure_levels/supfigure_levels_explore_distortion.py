import apcmodels.simulation as si
import apcmodels.anf as anf
import apcmodels.decode as dc
from apcmodels.util import save_to_csv
import numpy as np
import os, sys
sys.path.append(os.getcwd())
from util.functions import ISOToneGuest2021_exp1a, adjust_level
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Qt5Agg')
plt.ion()
from scipy.signal import butter, lfilter


def distortion_plot(F0=280, model=anf.AuditoryNerveZilany2014, model_name='Zilany2014', fs=300e3, lpf=True):
    # Define stimulus parameters
    levels = [0, 30, 60, 90]  # dB SPL, per component
    dur = 0.10  # seconds
    dur_ramp = 0.01  # seconds

    # Define model parameters
    cf_low = 5*F0
    cf_high = 11*F0
    n_cf = 40
    n_fiber_per_chan = round(((np.log10(11/5)/3)*18000)/n_cf)  # assume ~18k HSR fibers from 0.2 to 20 kHz

    # Encode parameters
    params = si.Parameters(dur=dur, dur_ramp=dur_ramp, fs=fs, n_cf=n_cf, delta_theta=[0.001], API=np.zeros(1),
                           n_fiber_per_chan=n_fiber_per_chan, model_name=model_name, F0=F0, cf_low=cf_low, cf_high=cf_high)
    params.wiggle('level', levels)                             # wiggle levels

    # Adjust levels to be in dB re: threshold
    params.flatten()
    for ele in params:
        ele['nominal_level'] = ele['level']                                 # encode nominal level (dB re: threshold)
        ele['level'] = adjust_level(ele['F0']*8, ele['level'], model_name)  # encode actual level (dB SPL)
    
    # Synthesize stimuli
    synth = ISOToneGuest2021_exp1a()
    stimuli = synth.synthesize_sequence(params)
    params.add_inputs(stimuli)

    # Construct simulation and run
    sim = model()
    results = sim.run(params)

    # Construct lowpass filter
    [b, a] = butter(2, 50000/(fs/2), 'low')
    if lpf:
        results = [lfilter(b, a, result, axis=1) for result in results]

    # Calculate FFT of each response along time axis
    spectra = [20*np.log10(np.abs(np.fft.fft(result, axis=1))) for result in results]

    # Construct plot, (1, n_level) subplots
    fig, axs = plt.subplots(len(levels), 1)
    f = np.linspace(0, fs, spectra[1].shape[1])
    for ii in range(len(spectra)):
        axs[ii].pcolormesh(f, np.arange(0, 40), spectra[ii])
        axs[ii].set_xlim((0, 150000))
        for harm in np.arange(1, 15):
            axs[ii].plot([harm*F0, harm*F0], [0, 40], 'r--')

# Notes:
# - I switched the sampling rate for the synapse stage up to 100 kHz in the local install of
# apcmodels... This seems to have eliminated tell-tale signs of aliasing in the response

distortion_plot(lpf=False)
distortion_plot(lpf=True)


