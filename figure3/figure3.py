"""
This script implements Figure 3 of Guest and Oxenham (2021).
"""

import apcmodels.simulation as si
import apcmodels.anf as anf
import numpy as np
from util.functions import ISOToneGuest2021, GEOMToneGuest2021
import matplotlib.pyplot as plt
import os, sys
sys.path.append(os.getcwd())


def plot_ep(axis_main, F0, condition, level, level_noise, title, first, color, fs=int(100e3)):
    """
    Function to plot excitation pattern of a single tone from Experiment 1.

    Arguments:
        axis_main (axis): axis object on which to plot
        F0 (float): F0 of the tone in Hz
        condition (str): condition from which to plot tones, either 'ISO' or 'GEOM'
        level (float): nominal level of the tones per-component, in dB SPL
        level_noise (float): nominal level of the noise, in units of dB SPL in the ERB centered at 1 kHz
        title (str): title to add to plot
        first (bool): whether or not this is the first plot in the row
        color (ndarray, str): the color to use for the plot
        fs (int): sampling rate to run the simulation at
    """
    # Figure out, for this F0, how many components will be in the complex tones
    n_freq = len(np.arange(F0, 48000 / 2, F0))
    n_freq_masker = len(np.arange(F0 * 2 ** (1 / 12), 48000 / 2, F0 * 2 ** (1 / 12)))
    # Set up simulation and stimulus parameters
    params = si.Parameters(F0=F0, F0_masker=F0 * 2 ** (1 / 12), fs=fs,
                           level=lambda: np.random.uniform(level-3, level+3, n_freq),
                           level_masker=lambda: np.random.uniform(level-3, level+3, n_freq_masker),
                           ten=True, level_noise=level_noise,
                           cf_low=F0 * 4, cf_high=F0 * 12, n_cf=200, fiber_type='msr')
    params.repeat(10)
    # Synthesize stimuli
    if condition == 'ISO':
        stimulus = ISOToneGuest2021().synthesize_sequence(params)
    else:
        stimulus = GEOMToneGuest2021().synthesize_sequence(params)
    # Add stimulus and flatten
    params.add_inputs(stimulus)
    params.flatten_and_unnest()
    # Estimate responses
    sim = anf.AuditoryNerveZilany2014()
    resp = sim.run(params)
    # Calculate cfs
    cfs = 10**np.linspace(np.log10(F0*4), np.log10(F0*12), 200)
    # Calculate mean over time and standard deivation over means
    firing_rates = np.array([np.mean(r, axis=1) for r in resp])
    mean_response = np.mean(firing_rates, axis=0)
    sd_response = np.std(firing_rates, axis=0)
    # Plot
    for harmonic in [6, 7, 8, 9, 10]:
        axis_main.plot([harmonic, harmonic], [0, 250], color='gray', linestyle='dashed', label='_nolegend_')
    axis_main.plot(cfs/F0, mean_response, color=color)
    axis_main.fill_between(cfs/F0, mean_response - sd_response, mean_response + sd_response,
                           color=color, alpha=0.2, label='_nolegend_')
    axis_main.set_xscale('log')
    axis_main.set_ylabel('Firing rate (sp/s)')
    if first:
        axis_main.get_xaxis().set_visible(False)
    if not first:
        axis_main.set_xlabel('CF (Harmonic Number)')
    axis_main.set_ylim((50, 250))
    axis_main.set_title(title)
    axis_main.set_xlim((5, 11))
    axis_main.set_xticks([5, 6, 7, 8, 9, 10, 11])
    axis_main.set_xticklabels(['', 6, 7, 8, 9, 10, ''])

# Plot excitation patterns
f, (a0, a1) = plt.subplots(2, 1, figsize=(4.25, 4))
plot_ep(a0, 280, 'ISO', 50, 40, '280 Hz', True, '#fc8d62')
plot_ep(a0, 280, 'GEOM', 50, 40, '280 Hz', True, '#8da0cb')
plot_ep(a1, 1400, 'ISO', 50, 40, '1400 Hz', False, '#fc8d62')
plot_ep(a1, 1400, 'GEOM', 50, 40, '1400 Hz', False, '#8da0cb')
plt.legend(['ISO', 'GEOM'], loc=3, framealpha=1)
plt.tight_layout()
# Save to disk
plt.savefig('plots/fig3.png', bbox_inches='tight')
