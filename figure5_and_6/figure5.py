"""
This script implements the simulations described in Figure 5a of Guest and Oxenham (2021).
"""

import apcmodels.simulation as si
import apcmodels.anf as anf
import numpy as np
import itertools
import os, sys
sys.path.append(os.getcwd())
from util.functions import DBLToneGuest2021
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import util as cfg
from util.functions import adjust_level


# Function to calculate autocorrelations
def calculate_autocorrelation(f0, neural_harm_nums, tmr, n_repeat=20, level_noise=37):
   # Setup params [simulation]
    params = si.Parameters(fs=int(100e3), fs_synapse=20e3, n_cf=1, anf_num=(1, 0, 0),                        # model
                           F0=f0, F0_masker_1=f0*2**(-5.5/12), F0_masker_2=f0*2**(6/12),                     # F0s
                           level=50, level_masker_1=50-tmr, level_masker_2=50-tmr,                           # levels
                           phase=0, phase_masker_1=0, phase_masker_2=0,                                      # phases
                           ten=True, level_noise=level_noise)                                                # noise
    params.wiggle_parallel(['cf_low', 'cf_high'], [neural_harm_nums*f0, neural_harm_nums*f0])

    # Synthesize stimuli and add to params
    stim = DBLToneGuest2021()
    params.add_inputs(stim.synthesize_sequence(params))

    # Adjust levels
    for ele in params:
        ele['level'] = adjust_level(ele['F0']*8, ele['level'], 'Zilany2014')  # encode adjusted level (dB SPL)

    # Encode repeats
    params.repeat(n_repeat)

    # Select model and run
    sim = anf.AuditoryNerveZilany2014()
    results = sim.run(params, runfunc=lambda x: [sim.simulate(ele) for ele in x])

    # Loop through CFs and repeats and transform firing rate into autocorrelation
    output = []
    for ele in results:
        output.append([np.fft.ifft(np.abs(np.fft.fft(np.squeeze(repeat)))**2) for repeat in ele])
    output = [np.mean(np.array(ele), axis=0) for ele in output]
    output = np.array(output)

    # Calculate lags
    lags = np.linspace(0, output.shape[1]/100e3, output.shape[1])

    # Get rid of extraneous samples
    keep = lags < (1/f0 * 20)
    lags = lags[keep]
    output = output[:, keep]

    # Return
    return lags, neural_harm_nums, np.abs(output)


# Run simulations and save to disk
tmrs = [0, 5, 10, 100] # set TMR for each simulation
f0s = [280, 1400]      # set F0s for each simulation

for tmr in tmrs:
    for f0 in f0s:
        lags, neural_harm_nums, histograms = calculate_autocorrelation(f0, np.linspace(4, 12, 100), tmr, n_repeat=10)
        np.save('figure5/autocorr_tmr_' + str(tmr) + '_' + str(f0) + '.npy', histograms)
        np.save('figure5/neural_harm_nums_tmr_' + str(tmr) + '_' + str(f0) + '.npy', neural_harm_nums)
        np.save('figure5/lags_tmr_' + str(tmr) + '_' + str(f0) + '.npy', lags)


# Functions to plot autocorrelation
def plot_autocorrelation(f0, xlims=None, cmap='viridis', figsize=None, tmrs=None, off_cf_only=False):
    """ Plots outputs of calculate_autocorrelation saved to disk

    Args:
        f0 (float): F0 of the stimulus in Hz
        xlims (None, tuple, list): either None (in which case default xlims are used) or a tuple/list containing
            xlims in units of periods of the F0

    Returns:
        fig (Figure): figure object generated, useful for adding to plot or adjusting settings post hoc
        ax (Axes): axes object generated, useful for adding to plot or adjusting settings post hoc
    """
    # Next, we handle the tmrs input argument
    if tmrs == None:
        tmrs = [0, 5, 10, 100]
    sacf_colors = {0: '#1f77b4', 5: '#ff7f0e', 10: '#2ca02c', 100: '#d62728'}
    sacf_colors = [sacf_colors[tmr] for tmr in tmrs]

    # First, we construct the figure 
    if figsize == None:
        fig, ax = plt.subplots(len(tmrs)+1, 1, figsize=(7, 2.25*(len(tmrs)+1)), sharex='all')
    else:
        fig, ax = plt.subplots(len(tmrs)+1, 1, figsize=figsize, sharex='all')
    
    # Now, we load in the simulations we want and plot each as we go
    for idx_tmr, tmr in enumerate(tmrs):
        # Load in histograms, neural harmonic numbres, and lags
        autocorr = np.load('figure5/autocorr_tmr_' + str(tmr) + '_' + str(f0) + '.npy')
        neural_harm_nums = np.load('figure5/neural_harm_nums_tmr_' + str(tmr) + '_' + str(f0) + '.npy')
        lags = np.load('figure5/lags_tmr_' + str(tmr) + '_' + str(f0) + '.npy')
        # If we want off_cf_only, we kill all channels outside of chosen
        if off_cf_only:
            chosen = np.sum(np.array([np.logical_and(0.3 < np.abs(neural_harm_nums - harm), np.abs(neural_harm_nums - harm) < 1) for harm in [6, 7, 8, 9, 10]]), axis=0) > 1
            autocorr[np.logical_not(chosen), :] = np.NaN
        # Plot surface color plot using pcolormesh
        ax[idx_tmr].pcolormesh(lags*f0, neural_harm_nums, autocorr, cmap=cmap)
        # Set labels
        if xlims == None:
            ax[idx_tmr].set_xlim((0, 8))
        else:
            ax[idx_tmr].set_xlim(xlims)
        # Hide x-axis 
        ax[idx_tmr].get_xaxis().set_visible(False)
        # Set y-ticks
        ax[idx_tmr].set_yticks([6, 7, 8, 9, 10])
        # Calculate average across simulated CFs and plot on bottom row
        means = np.nanmean(autocorr, axis=0)
        ax[len(tmrs)].plot(lags*f0, means/np.max(means), color=sacf_colors[idx_tmr])
        # Set labels
        if xlims == None:
            ax[len(tmrs)].set_xlim((0, 8))
        else:
            ax[len(tmrs)].set_xlim(xlims)
        ax[len(tmrs)].get_yaxis().set_visible(False)
        ax[len(tmrs)].set_xlabel('Lag (periods of F0)')
    # Save to disk
    plt.tight_layout()
    # Return figures
    return fig, ax


# Create plots
# Define colorscheme for different F0s
color_lower = '#ff3b00'
color_middle = '#ffb200'
color_upper = '#ff00ac'
linewidth = 2
linestyle = 'dashed'
cmap = 'viridis'

# Create main plot (top part)
for fig_subnum, F0 in zip(['a', 'b'], [280, 1400]):
    fig, ax = plot_autocorrelation(F0, cmap=cmap, figsize=(7*0.7, 12*0.7))
    # Adjust t-ticks
    for idx_tmr in range(4):
        ax[idx_tmr].set_yticks([6, 7, 8, 9, 10])
    plt.savefig('plots/fig5' + fig_subnum + '1.png', dpi=300)
plt.close('all')

# Create zoomed-in versions of plots
# #1: Low frequencies, from 3.5 to 4.5, only 0/5/10 dB
fig, ax = plot_autocorrelation(280, (3.5, 4.5), cmap=cmap, figsize=(5*0.6, 10*0.6), tmrs=[0, 5, 10])
for idx_tmr in range(3):
    # Highlight peaks at F0 multiples
    for mult in np.arange(1, 11):
        # Plot at F0 period
        ax[idx_tmr].plot([mult, mult], [4, 12], color=color_middle, linewidth=linewidth, linestyle=linestyle)
        # Plot at masker period #1
        ax[idx_tmr].plot([mult/2**(-5.5/12), mult/2**(-5.5/12)], [4, 12], color=color_lower, linewidth=linewidth, linestyle=linestyle)
        # Plot at masker period #2
        ax[idx_tmr].plot([mult/2**(6/12), mult/2**(6/12)], [4, 12], color=color_upper, linewidth=linewidth, linestyle=linestyle)
plt.savefig(os.path.join('plots', 'fig5' + 'zoom1.png'), dpi=300)

# #2: High frequencies, from 1.5 to 2.5, only 5/10/quiet dB
fig, ax = plot_autocorrelation(1400, (1.5, 2.5), cmap=cmap, figsize=(5*0.6, 10*0.6), tmrs=[5, 10, 100])
for idx_tmr in range(3):
    # Highlight peaks at F0 multiples
    for mult in np.arange(1, 11):
        # Plot at F0 period
        ax[idx_tmr].plot([mult, mult], [4, 12], color=color_middle, linewidth=linewidth, linestyle=linestyle)
        if idx_tmr < 2:
            # Plot at masker period #1
            ax[idx_tmr].plot([mult/2**(-5.5/12), mult/2**(-5.5/12)], [4, 12], color=color_lower, linewidth=linewidth, linestyle=linestyle)
            # Plot at masker period #2
            ax[idx_tmr].plot([mult/2**(6/12), mult/2**(6/12)], [4, 12], color=color_upper, linewidth=linewidth, linestyle=linestyle)
plt.savefig(os.path.join('plots', 'fig5' + 'zoom2.png'), dpi=300)


# #3: Same as #2, but off-cf fibers only in sACF
fig, ax = plot_autocorrelation(1400, (1.5, 2.5), cmap=cmap, figsize=(5*0.6, 10*0.6), tmrs=[5, 10, 100], off_cf_only=True)
for idx_tmr in range(3):
    # Highlight peaks at F0 multiples
    for mult in np.arange(1, 11):
        # Plot at F0 period
        ax[idx_tmr].plot([mult, mult], [4, 12], color=color_middle, linewidth=linewidth, linestyle=linestyle)
        # Plot at masker period #1
        if idx_tmr < 2:
            ax[idx_tmr].plot([mult/2**(-5.5/12), mult/2**(-5.5/12)], [4, 12], color=color_lower, linewidth=linewidth, linestyle=linestyle)
            # Plot at masker period #2
            ax[idx_tmr].plot([mult/2**(6/12), mult/2**(6/12)], [4, 12], color=color_upper, linewidth=linewidth, linestyle=linestyle)
plt.savefig(os.path.join('plots', 'fig5' + 'zoom3.png'), dpi=300)