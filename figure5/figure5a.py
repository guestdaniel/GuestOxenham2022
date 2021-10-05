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


def calculate_ISI_histogram(f0, neural_harm_nums, tmr, n_repeat=20):
    """ Synthesizes complex tone mixture and simulates interspike interval histogram response.

    Args:
        f0 (float): cf of the complex tone
        neural_harm_nums (np.ndarray): array of neural harmonic numbers (CF/F0) to test. This is performed by keeping
            the F0 of the tone fixed while manipulating the range of CFs tested.
        tmr (float): Target-to-masker ratio of the target complex tone vs the masker complex tone
        n_repeat (int): number of repeats to run for each CF/F0 to test

    Returns:
        histograms (list): list of histograms counts for ISI histogram for each neural harmonic number
    """
    # Setup params [simulation]
    params = si.Parameters(fs=int(200e3), n_cf=1, anf_num=(1, 0, 0),                                          # model
                           F0=280, F0_masker_1=280*2**(-5.5/12), F0_masker_2=280*2**(6/12),                   # F0s
                           level=50, level_masker_1=50-tmr, level_masker_2=50-tmr,                            # levels
                           phase=0, phase_masker_1=0, phase_masker_2=0,                                       # phases
                           ten=True, level_noise=40)                                                          # noise
    params.wiggle_parallel(['cf_low', 'cf_high'], [neural_harm_nums*f0, neural_harm_nums*f0])

    # Synthesize stimuli and add to params
    stim = DBLToneGuest2021()
    params.add_inputs(stim.synthesize_sequence(params))

    # Encode repeats
    params.repeat(n_repeat)

    # Select model and run
    sim = anf.AuditoryNerveZilany2014Spikes()
    results = sim.run(params, runfunc=lambda x: [sim.simulate(ele) for ele in x])

    # Extract spike times from results
    histograms = []  # top-level list that will store the histogram for all neural harm numbers
    for idx_result, result in enumerate(results):
        isis_this_channel = []  # list will store the isis for this neural harm number
        for repeat in result:
            # Extract spikes from data frame
            spike_times = repeat.loc[0]['spikes']
            # Calculate ISIs
            isi = [np.abs(spike_times[idx[1]] - spike_times[idx[0]]) for idx in
                         itertools.combinations(list(range(len(spike_times))), r=2)]
            isis_this_channel.append(isi)
        # Combine ISIs into one list/array
        ISIs = list(itertools.chain(*isis_this_channel))
        # Transform ISIs from seconds into units of F0 period
        ISIs_period = np.array(ISIs) / (1 / f0)
        # Calculate histogram and append to results
        hist_low, edges_low = np.histogram(ISIs_period, bins=2200, range=(0, 20))
        histograms.append(hist_low)

    # Return
    return histograms


def plot_ISI_histogram(f0, xlims=None, cmap='viridis'):
    """ Plots outputs of calculate_ISI_histogram saved to disk

    Args:
        f0 (float): F0 of the stimulus in Hz
        xlims (None, tuple, list): either None (in which case default xlims are used) or a tuple/list containing
            xlims in units of periods of the F0

    Returns:
        fig (Figure): figure object generated, useful for adding to plot or adjusting settings post hoc
        ax (Axes): axes object generated, useful for adding to plot or adjusting settings post hoc
    """
    # Create main plot (part a)
    fig, ax = plt.subplots(4, 1, figsize=(7, 9), sharex='all')
    for idx_tmr, tmr in enumerate(tmrs):
        # Load in histograms and calculate edges
        histograms = np.load('figure5/isi_histograms_tmr_' + str(tmr) + '_' + str(f0) + '.npy')
        neural_harm_nums = np.load('figure5/neural_harm_nums_tmr_' + str(tmr) + '_' + str(f0) + '.npy')
        edges = np.linspace(start=0, stop=20, num=2200)
        # Plot pcolormesh
        ax[idx_tmr].pcolormesh(edges, neural_harm_nums, histograms, cmap=cmap)
        # Set labels
        if xlims == None:
            ax[idx_tmr].set_xlim((0, 8))
        else:
            ax[idx_tmr].set_xlim(xlims)
        # Hide x-axis
        if idx_tmr > 0:
            ax[idx_tmr].get_xaxis().set_visible(False)
        ax[idx_tmr].xaxis.set_ticks_position('top')
    # Save to disk
    plt.tight_layout()
    # Return figures
    return fig, ax


# Run simulations
tmrs = [0, 4, 8, 12]  # set TMR for each simulation
f0s = [280, 1400]     # set F0s for each simulation

for tmr in tmrs:
    for f0 in f0s:
        # Calculate ISI histograms for one Larsen and Delgutte (2008) stimulus and save to disk
        histograms = calculate_ISI_histogram(f0, np.linspace(4, 12, 40), tmr=tmr, n_repeat=100)
        np.save('figure5/isi_histograms_tmr_' + str(tmr) + '_' + str(f0) + '.npy', histograms)
        np.save('figure5/neural_harm_nums_tmr_' + str(tmr) + '_' + str(f0) + '.npy', np.linspace(4, 12, 40))

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