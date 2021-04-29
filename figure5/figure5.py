"""
This script implements the simulations described in Figure 5 of Guest and Oxenham (2021).
"""
import apcmodels.simulation as si
import apcmodels.anf as anf
import numpy as np
from util.functions import DBLToneGuest2021
import matplotlib.pyplot as plt
import os, sys
sys.path.append(os.getcwd())
from functools import partial


def prep_ep(F0, level, level_maskers, level_noise, fs=int(100e3)):
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
    masker_F0_1 = F0 * 2 ** (-6.25/12)
    masker_F0_2 = F0 * 2 ** (5.25/12)

    # Figure out, for this F0, how many components will be in the complex tones
    n_freq_masker_1 = len(np.arange(masker_F0_1, 48000 / 2, masker_F0_1))
    n_freq_masker_2 = len(np.arange(masker_F0_2, 48000 / 2, masker_F0_2))

    # Set up simulation and stimulus parameters
    params = si.Parameters(F0_masker_1=masker_F0_1, F0_masker_2=masker_F0_2,                                     # F0s
                           level_masker_1=lambda: np.random.uniform(level_maskers-3, level+3, n_freq_masker_1),  # levels
                           level_masker_2=lambda: np.random.uniform(level_maskers-3, level+3, n_freq_masker_2),  # levels
                           ten=True, level_noise=level_noise,                                                    # noise
                           cf_low=F0 * 4, cf_high=F0 * 12, n_cf=200, fiber_type='msr',                           # anf
                           fs=fs)

    return params


def handle_labels_and_axes(axis_main, first, title):
    """ Helper function for plotting excitation patterns

    Args:
        axis_main (axis): axis object on which to plot
        title (str): title to add to plot
        first (bool): whether or not this is the first plot in the row
    """
    axis_main.set_xscale('log')
    axis_main.set_ylabel('Firing rate (sp/s)')
    if first:
        axis_main.get_xaxis().set_visible(False)
    if not first:
        axis_main.set_xlabel('CF (Harmonic Number)')
    axis_main.set_title(title)
    axis_main.set_xlim((5, 11))
    axis_main.set_xticks([5, 6, 7, 8, 9, 10, 11])
    axis_main.set_xticklabels(['', 6, 7, 8, 9, 10, ''])


def plot_ep(axis_main, F0, level, level_maskers, level_noise, title, first, color, fs=int(100e3)):
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
    """
    # Get params
    params = prep_ep(F0, level, level_maskers, level_noise, fs)
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
    # Plot
    for harmonic in [6, 7, 8, 9, 10]:
        axis_main.plot([harmonic, harmonic], [0, 250], color='gray', linestyle='dashed', label='_nolegend_')
    axis_main.plot(cfs/F0, mean_response, color=color)
    axis_main.fill_between(cfs/F0, mean_response - sd_response, mean_response + sd_response,
                           color=color, alpha=0.2, label='_nolegend_')
    # Handle labels and axes
    handle_labels_and_axes(axis_main, first, title)
    axis_main.set_ylim((50, 250))


def plot_ep_difference(axis_main, F0, inc_size, level, level_maskers, level_noise, title, first, color, fs=int(100e3)):
    """ Function to plot excitation pattern of DBL tones from Experiment 2

    Args:
        axis_main (axis): axis object on which to plot
        F0 (float): F0 of the tone in Hz. The two maskers tones are set to +/- 6.25 ST
        inc_size (float): size of the increment in F0 between the first interval and second interval, in ST
        level (float): level per-component in dB SPL of the targets before filtering
        level_maskers (float); level per-component in dB SPL of the maskers before filtering
        level_noise (float): level of the TEN in the ERB centered on 1 kHz
        title (str): title to add to plot
        first (bool): whether or not this is the first plot in the row
        color (ndarray, str): the color to use for the plot
        fs (int): sampling rate, in Hz
    """
    # Get params
    params = prep_ep(F0, level, level_maskers, level_noise, fs)
    params.wiggle_parallel(['F0', 'level'], [[F0, F0*2**(inc_size/12)],
                                             [lambda: np.random.uniform(level-3, level+3, len(np.arange(F0, 48000 / 2, F0))),
                                              lambda: np.random.uniform(level-3, level+3, len(np.arange(F0, 48000 / 2, F0*2**(inc_size/12))))]])
    params.repeat(10)
    # Synthesize stimuli
    stimulus = DBLToneGuest2021().synthesize_sequence(params)
    # Add stimulus
    params.add_inputs(stimulus)
    # Define special runfunc needed for this simulation
    def run_elements(sim, paramlist):
        return [sim.simulate(x) for x in paramlist]
    # Estimate responses
    sim = anf.AuditoryNerveZilany2014()
    resp = sim.run(params, runfunc=partial(run_elements, sim))
    # Calculate cfs
    cfs = 10**np.linspace(np.log10(F0*4), np.log10(F0*12), 200)
    # For each simulation pair, calculate the difference between sims from reference and incremented F0
    differences = np.array([np.mean(x[1], axis=1) - np.mean(x[0], axis=1) for x in zip(resp[0], resp[1])])
    # Calculate mean over time and standard deivation over means
    mean_response = np.mean(differences, axis=0)
    sd_response = np.std(differences, axis=0)
    # Plot arrow markers at harmonics of the F0 and harmonics of the target
    for harmonic in [6, 7, 8, 9, 10]:
        # Arrows
        axis_main.arrow(harmonic - 0.11, mean_response[np.argmin(np.abs(cfs/F0 - harmonic))] - 35, 0.1, 25, color='black', width=0.025,
                  head_width=0.1, head_length=6, length_includes_head=True)
        axis_main.arrow(harmonic * 2 ** (inc_size/12) + 0.11, mean_response[np.argmin(np.abs(cfs/F0 - harmonic * 2 ** (inc_size/12)))] + 35, -0.1, -25, color='red', width=0.025,
                  head_width=0.1, head_length=6, length_includes_head=True)
    # Plot simulation data
    axis_main.plot(cfs/F0, mean_response, color=color)
    axis_main.fill_between(cfs/F0, mean_response - sd_response, mean_response + sd_response,
                           color=color, alpha=0.2, label='_nolegend_')
    # Handle labels and axes
    handle_labels_and_axes(axis_main, first, title)
    axis_main.set_ylim((-100, 100))


# Plot excitation patterns [MAIN]
f, (a0, a1) = plt.subplots(2, 1, figsize=(4.25, 4))
plot_ep(a0, 280, 50, 50, 40, '280 Hz', True, '#7fc97f')
plot_ep(a1, 1400, 50, 50, 40, '1400 Hz', False, '#a28ac1')
plt.tight_layout()
plt.savefig('plots/fig5a.png', bbox_inches='tight')

# Plot excitation patterns [MAIN]
f, (a0, a1) = plt.subplots(2, 1, figsize=(4.25, 4))
plot_ep(a0, 280, 50, 44, 40, '280 Hz', True, '#7fc97f')
plot_ep(a1, 1400, 50, 39, 40, '1400 Hz', False, '#a28ac1')
plt.tight_layout()
plt.savefig('plots/fig5b.png', bbox_inches='tight')

# Plot excitation patterns [MAIN]
f, (a0, a1) = plt.subplots(2, 1, figsize=(4.25, 4))
plot_ep_difference(a0, 280, 1, 50, 44, 40, '280 Hz', True, '#7fc97f')
plot_ep_difference(a1, 1400, 1, 50, 39, 40, '1400 Hz', False, '#a28ac1')
plt.tight_layout()
plt.savefig('plots/fig5c.png', bbox_inches='tight')