"""
This script implements the simulations described in Figure 5b of Guest and Oxenham (2021).
"""
import apcmodels.simulation as si
import apcmodels.anf as anf
import numpy as np
from util.functions import DBLToneGuest2021
import matplotlib.pyplot as plt
import os, sys
sys.path.append(os.getcwd())
from functools import partial


def prep_ep(F0, level, level_maskers, level_noise, fs=int(100e3), fiber_type='msr'):
    """ Helper function for plotting excitation patterns

    Args:
        F0 (float): F0 of the tone in Hz. The two maskers tones are set to +/- 6.25 ST
        level (float): level per-component in dB SPL of the targets before filtering
        level_maskers (float); level per-component in dB SPL of the maskers before filtering
        level_noise (float): level of the TEN in the ERB centered on 1 kHz
        fs (int): sampling rate, in Hz

    Returns:
        params (Parameters): Parameters object with everything except for the F0 encoded in the simulation
    """
    # Set masker F0s
    masker_F0_1 = F0 * 2 ** (-5.5/12)
    masker_F0_2 = F0 * 2 ** (6/12)

    # Figure out, for this F0, how many components will be in the complex tones
    n_freq_masker_1 = len(np.arange(masker_F0_1, 48000 / 2, masker_F0_1))
    n_freq_masker_2 = len(np.arange(masker_F0_2, 48000 / 2, masker_F0_2))

    # Set up simulation and stimulus parameters
    params = si.Parameters(F0_masker_1=masker_F0_1, F0_masker_2=masker_F0_2,
                           level_masker_1=lambda: np.random.uniform(level_maskers-3, level+3, n_freq_masker_1),  # levels
                           level_masker_2=lambda: np.random.uniform(level_maskers-3, level+3, n_freq_masker_2),  # levels
                           ten=True, level_noise=level_noise,                                                    # noise
                           cf_low=F0 * 4, cf_high=F0 * 12, n_cf=200, fiber_type=fiber_type,                      # anf
                           fs=fs)

    return params


def handle_labels_and_axes(axis_main, first, title, yaxis_side='left'):
    """ Helper function for plotting excitation patterns

    Args:
        axis_main (axis): axis object on which to plot
        title (str): title to add to plot
        first (bool): whether or not this is the first plot in the row
    """
    axis_main.set_xscale('log')
    #axis_main.set_ylabel('Firing rate (sp/s)')
    if first:
        axis_main.get_xaxis().set_visible(False)
    #if not first:
        #axis_main.set_xlabel('CF (Harmonic Number)')
    #axis_main.set_title(title)
    axis_main.set_xlim((4, 12))
    axis_main.set_xticks([4, 5, 6, 7, 8, 9, 10, 11, 12])
    axis_main.set_xticklabels(['', '', 6, 7, 8, 9, 10, '', ''])
    if yaxis_side == 'right':
        axis_main.yaxis.set_ticks_position('right')



def plot_ep(axis_main, F0, level, level_maskers, level_noise, title, first, color, fs=int(100e3), fiber_type='msr', yaxis_side='left'):
    """ Function to plot excitation pattern of DBL tones from Experiment 2

    Args:
        axis_main (axis): axis object on which to plot
        F0 (float): F0 of the tone in Hz. The two maskers tones are set to +/- 6.25 ST
        level (float): level per-component in dB SPL of the targets before filtering
        level_maskers (float); level per-component in dB SPL of the maskers before filtering
        level_noise (float): level of the TEN in the ERB centered on 1 kHz
        title (str): title to add to plot
        first (bool): whether or not this is the first plot in the row
        color (ndarray, str): the color to use for the plot
        fs (int): sampling rate, in Hz
        fiber_type (str): which fiber type ('high', 'medium', or 'low') to use
        yaxis_side (str): where to plot the yticks and yticklabels ('left' or 'right')
    """
    # Get params
    params = prep_ep(F0, level, level_maskers, level_noise, fs, fiber_type)
    params.append(['F0', 'level'], [F0, lambda: np.random.uniform(level-3, level+3, len(np.arange(F0, 48000 / 2, F0)))])
    params.repeat(10)
    # Synthesize stimuli
    stimulus = DBLToneGuest2021().synthesize_sequence(params)
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
    # Plot masker components
    for harmonic in [2, 3, 4, 5, 6, 7, 8, 8, 10, 11, 12]:
        idx_min = np.argmin(np.abs(cfs/F0 - harmonic*2**(6/12)))
        axis_main.plot([harmonic*2**(6/12), harmonic*2**(6/12)], [0, mean_response[idx_min]], color='cyan',
                       linestyle='dashed', label='_nolegend_')
        idx_min = np.argmin(np.abs(cfs/F0 - harmonic*2**(-5.5/12)))
        axis_main.plot([harmonic*2**(-5.5/12), harmonic*2**(-5.5/12)], [0, mean_response[idx_min]], color=[0, 1, 0],
                       linestyle='dashed', label='_nolegend_')
    # Plot target components
    for harmonic in [6, 7, 8, 9, 10]:
        idx_min = np.argmin(np.abs(cfs/F0 - harmonic))
        axis_main.plot([harmonic, harmonic], [0, mean_response[idx_min]], color='blue', linestyle='dashed',
                       label='_nolegend_', linewidth=2)
    axis_main.plot(cfs/F0, mean_response, color=color)
    axis_main.fill_between(cfs/F0, mean_response - sd_response, mean_response + sd_response,
                           color=color, alpha=0.2, label='_nolegend_')
    # Handle labels and axes
    handle_labels_and_axes(axis_main, first, title, yaxis_side)
    if fiber_type == 'hsr':
        axis_main.set_ylim((0, 375))
    elif fiber_type == 'msr':
        axis_main.set_ylim((0, 250))
    else:
        axis_main.set_ylim((0, 50))


# Plot excitation patterns for HSR and LSR fibers at various TMRs
for F0, label, side in zip([280, 1400], ['c', 'd'], ['left', 'right']):
    f, axs = plt.subplots(4, 2, figsize=(6, 6), sharex='all')
    tmrs = [0, 4, 8, 12]
    fiber_types = ['hsr', 'lsr']
    for idx_tmr in range(4):
        for idx_fiber_type in range(2):
            plot_ep(axs[idx_tmr, idx_fiber_type], F0, 50, 50 - tmrs[idx_tmr], 40, '1400 Hz', False, '#a28ac1',
                    fiber_type=fiber_types[idx_fiber_type], yaxis_side=side)
    plt.savefig('plots/fig5' + label + '1.png', bbox_inches='tight')