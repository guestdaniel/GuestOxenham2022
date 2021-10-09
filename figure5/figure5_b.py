"""
This script implements the simulations described in Figure 5b of Guest and Oxenham (2021).
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


def summarize_isi_histograms(f0, tolerance=0.05):
    """ Plots outputs of calculate_ISI_histogram saved to disk

    Args:
        f0 (float): F0 of the stimulus in Hz
        xlims (None, tuple, list): either None (in which case default xlims are used) or a tuple/list containing
            xlims in units of periods of the F0

    Returns:
        fig (Figure): figure object generated, useful for adding to plot or adjusting settings post hoc
        ax (Axes): axes object generated, useful for adding to plot or adjusting settings post hoc
    """
    # Load data in
    histograms = []
    neural_harm_nums = np.load('figure5/neural_harm_nums_tmr_' + str(0) + '_' + str(f0) + '.npy')
    edges = np.linspace(start=0, stop=20, num=2200)
    for idx_tmr, tmr in enumerate(tmrs):
        # Load in histograms and calculate edges
        histograms.append(np.load('figure5/isi_histograms_tmr_' + str(tmr) + '_' + str(f0) + '.npy'))
    # Compile data across reps 
    # To summarize the isi histograms, we can compute some sort of signal to noise ratio --- maybe number of spikes observed near multiples of F0 period as a percentage of total number of spikes observed?
    def calculate_metric(histogram):
        harmonics = np.arange(1, 21)
        total_spikes = sum(sum(histogram))
        spikes_within_template = 0
        for harm in harmonics:
            spikes_within_template += sum(sum(histogram[:, np.logical_and(edges > (harm-tolerance), edges < (harm+tolerance))]))
        return spikes_within_template/total_spikes
    # Apply metric to all histograms
    results = [calculate_metric(x) for x in histograms]
    # Visualize
    plt.plot(tmrs, results)
    plt.ylim((0, 1))
    plt.tight_layout()

plt.figure()
summarize_isi_histograms(280, tolerance=0.01)
summarize_isi_histograms(1400, tolerance=0.01)

# Run simulations
tmrs = [0, 4, 8, 12]  # set TMR for each simulation
f0s = [280, 1400]     # set F0s for each simulation

# Create plots
# Define colorscheme for different F0s
color_lower = '#ff3b00'
color_middle = '#ffb200'
color_upper = '#ff00ac'
linewidth = 4
linestyle = 'dashed'
cmap = 'viridis'

# Create main plot (zoomed-out histograms surfaces)
for fig_subnum, F0 in zip(['a', 'b'], [280, 1400]):
    plot_ISI_histogram(F0, cmap=cmap)
    plt.savefig('plots/fig5' + fig_subnum + '1.png')

    # Create subplot #1 (zoomed-in histogram surface, focused on 3.5 to 4.5 periods of F0)
    fig, ax = plot_ISI_histogram(F0, (3.5, 4.5), cmap=cmap)
    for idx_tmr in range(4):
        # Highlight peaks at F0 multiples
        ax[idx_tmr].plot([4, 4], [4, 12], color=color_middle, linewidth=linewidth, linestyle=linestyle)
        ax[idx_tmr].plot([3/2**(-5.5/12), 3/2**(-5.5/12)], [4, 12], color=color_lower, linewidth=linewidth, linestyle=linestyle)
        ax[idx_tmr].plot([5/2**(6/12), 5/2**(6/12)], [4, 12], color=color_upper, linewidth=linewidth, linestyle=linestyle)
        ax[idx_tmr].plot([6/2**(6/12), 6/2**(6/12)], [4, 12], color=color_upper, linewidth=linewidth, linestyle=linestyle)
    plt.savefig('plots/fig5' + fig_subnum + '2.png')

    # Create subplot #2 (zoomed-in histogram surface, focused on 0.5 to 1.5 periods of F0)
    fig, ax = plot_ISI_histogram(F0, (0.5, 1.5), cmap=cmap)
    for idx_tmr in range(4):
        # Highlight peaks at F0 multiples
        ax[idx_tmr].plot([1, 1], [4, 12], color=color_middle, linewidth=linewidth, linestyle=linestyle)
        ax[idx_tmr].plot([1/2**(-5.5/12), 1/2**(-5.5/12)], [4, 12], color=color_lower, linewidth=linewidth, linestyle=linestyle)
        ax[idx_tmr].plot([1/2**(6/12), 1/2**(6/12)], [4, 12], color=color_upper, linewidth=linewidth, linestyle=linestyle)
        ax[idx_tmr].plot([2/2**(6/12), 2/2**(6/12)], [4, 12], color=color_upper, linewidth=linewidth, linestyle=linestyle)
    plt.savefig('plots/fig5' + fig_subnum + '3.png')
plt.close('all')