import apcmodels.synthesis as sy
import apcmodels.simulation as si
import apcmodels.anf as anf
import apcmodels.decode as dc
from apcmodels.util import save_to_csv
import numpy as np
import os, sys
sys.path.append(os.getcwd())
from util.functions import ISOToneGuest2021_exp1a, adjust_level
import util as cfg
import matplotlib.pyplot as plt

def simulate_complex_tone_response(F0, F0_base, level, model=anf.AuditoryNerveZilany2014, 
                                model_name='Zilany2014', fs=200e3, 
                                fs_synapse=100e3, n_cf=40):
    # Define stimulus parameters
    dur = 0.10  # seconds
    dur_ramp = 0.01  # seconds

    # Define model parameters
    cf_low = 5*F0_base
    cf_high = 11*F0_base
    cfs = 10**np.linspace(np.log10(cf_low), np.log10(cf_high), n_cf)
    n_cf = n_cf 
    n_fiber_per_chan = round(((np.log10(11/5)/3)*18000)/n_cf)  # assume ~18k HSR fibers from 0.2 to 20 kHz

    # Encode parameters
    params = si.Parameters(dur=dur, dur_ramp=dur_ramp, fs=fs, delta_theta=[0.001], API=np.zeros(1),
                           n_fiber_per_chan=n_fiber_per_chan, model_name=model_name, F0=F0, level=level,
                           cfs=cfs, fs_synapse=fs_synapse)

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

    return cfs, results[0]


def simulate_pure_tone_response(freq, freq_base, level, model=anf.AuditoryNerveZilany2014, 
                                model_name='Zilany2014', fs=200e3, 
                                fs_synapse=100e3, n_cf=40):
    # Define stimulus parameters
    dur = 0.10  # seconds
    dur_ramp = 0.01  # seconds

    # Define model parameters
    cf_low = 0.25*freq_base
    cf_high = 1.75*freq_base
    cfs = 10**np.linspace(np.log10(cf_low), np.log10(cf_high), n_cf)
    n_fiber_per_chan = round(((np.log10(1.5/0.5)/3)*18000)/n_cf)  # assume ~18k HSR fibers from 0.2 to 20 kHz

    # Encode parameters
    params = si.Parameters(dur=dur, dur_ramp=dur_ramp, fs=fs, delta_theta=[0.001], API=np.zeros(1),
                           n_fiber_per_chan=n_fiber_per_chan, model_name=model_name, freq=freq, level=level,
                           cfs=cfs, fs_synapse=fs_synapse)

    # Adjust levels to be in dB re: threshold
    params.flatten()
    for ele in params:
        ele['nominal_level'] = ele['level']                                 # encode nominal level (dB re: threshold)
        ele['level'] = adjust_level(ele['freq'], ele['level'], model_name)  # encode actual level (dB SPL)

    # Synthesize stimuli
    synth = sy.PureTone()
    stimuli = synth.synthesize_sequence(params)
    params.add_inputs(stimuli)

    # Construct simulation and run
    sim = model()
    results = sim.run(params)
    return cfs, results[0]

# FDLs
freqs = 8 * 10**np.linspace(np.log10(280) - 0.2, np.log10(1400) + 0.1, 24)  # simulate 8th harmonic of F0s
levels = [20, 40, 60, 80]
fig, axs = plt.subplots(1, 3, figsize=(8, 2.5), sharey=True)
for idx, idx_freq in zip([0, 1, 2], [0, 10, 20]):
    for level in levels:
        cfs, x = simulate_pure_tone_response(freqs[idx_freq], freqs[idx_freq], level, n_cf=100, fs=100e3, fs_synapse=20e3)
        axs[idx].plot(cfs/freqs[idx_freq], np.mean(x, axis=1))
        axs[idx].set_title('Freq = ' + str(round(freqs[idx_freq])) + ' Hz')
    axs[idx].set_xlabel('CF (Hz)')
    axs[idx].plot([0.5, 0.5], [0, 300], color='gray', linestyle='dashed')
    axs[idx].plot([0.6, 0.6], [0, 300], color='teal', linestyle='dashed')
    axs[idx].plot([1.5, 1.5], [0, 300], color='gray', linestyle='dashed')
axs[0].set_ylabel('Firing rate (sp/s)')
plt.savefig(os.path.join('plots', 'supfigure_distortion_a.png'), dpi=300)

# F0DLs
F0s = 10**np.linspace(np.log10(280) - 0.2, np.log10(1400) + 0.1, 24)  # simulate 8th harmonic of F0s
levels = [20, 40, 60, 80]
fig, axs = plt.subplots(1, 3, figsize=(8, 2.5), sharey=True)
for idx, idx_F0 in zip([0, 1, 2], [0, 10, 20]):
    for level in levels:
        cfs, x = simulate_complex_tone_response(F0s[idx_F0], F0s[idx_F0], level, n_cf=100, fs=100e3, fs_synapse=20e3)
        axs[idx].plot(cfs/F0s[idx_F0], np.mean(x, axis=1))
        axs[idx].set_title('F0 = ' + str(round(F0s[idx_F0])) + ' Hz')
    axs[idx].set_xlabel('CF (Hz)')
    axs[idx].plot([5, 5], [0, 300], color='gray', linestyle='dashed')
    axs[idx].plot([11, 11], [0, 300], color='gray', linestyle='dashed')
axs[0].set_ylabel('Firing rate (sp/s)')
plt.savefig(os.path.join('plots', 'supfigure_distortion_b.png'), dpi=300)