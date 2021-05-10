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


def calculate_ISI_histogram(f0, neural_harm_nums, tmr):
    """ Synthesizes complex tone mixture and simulates interspike interval histogram response.

    Args:
        f0 (float): cf of the complex tone
        neural_harm_nums (np.ndarray): array of neural harmonic numbers (CF/F0) to test. This is performed by keeping
            the F0 of the tone fixed while manipulating the range of CFs tested.
        tmr (float): Target-to-masker ratio of the target complex tone vs the masker complex tone

    Returns:
        histograms (list): list of histograms counts for ISI histogram for each neural harmonic number
    """
    # Setup params [simulation]
    params = si.Parameters(fs=int(200e3), n_cf=1, anf_num=(0, 1, 0),                                          # model
                           F0=280, F0_masker_1=280*2**(-5.5/12), F0_masker_2=280*2**(6/12),                   # F0s
                           level=50, level_masker_1=50-tmr, level_masker_2=50-tmr,                            # levels
                           phase=0, phase_masker_1=0, phase_masker_2=0,                                       # phases
                           ten=True, level_noise=40)                                                          # noise
    params.wiggle_parallel(['cf_low', 'cf_high'], [neural_harm_nums*f0, neural_harm_nums*f0])

    # Synthesize stimuli and add to params
    stim = DBLToneGuest2021()
    params.add_inputs(stim.synthesize_sequence(params))

    # Encode repeats
    params.repeat(200)

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

for tmr in [0, 5, 10, 15]:
    # Calculate ISI histograms for one Larsen and Delgutte (2008) stimulus and save to disk
    histograms = calculate_ISI_histogram(280, np.linspace(4, 12, 40), tmr=tmr)
    np.save('figure5/isi_histograms_tmr_' + str(tmr) + '.npy', histograms)
    np.save('figure5/neural_harm_nums_tmr_' + str(tmr) + '.npy', np.linspace(4, 12, 40))

# Create main plot (part a)
fig, ax = plt.subplots(4, 1, figsize=(7, 9), sharex='all')
for idx_tmr, tmr in enumerate([0, 5, 10, 15]):
    # Load in histograms and calculate edges
    histograms = np.load('figure5/isi_histograms_tmr_' + str(tmr) + '.npy')
    neural_harm_nums = np.load('figure5/neural_harm_nums_tmr_' + str(tmr) + '.npy')
    edges = np.linspace(start=0, stop=20, num=2200)
    # Plot pcolormesh
    ax[idx_tmr].pcolormesh(edges, neural_harm_nums, histograms, cmap='plasma')
    # Set labels
    ax[idx_tmr].set_xlim((0, 8))
    # Hide x-axis
    if idx_tmr > 0:
        ax[idx_tmr].get_xaxis().set_visible(False)
    ax[idx_tmr].xaxis.set_ticks_position('top')
# Save to disk
plt.tight_layout()
plt.savefig('plots/fig5a1.png')

# Create subplot (part b)
fig, ax = plt.subplots(4, 1, figsize=(4, 6), sharex='all')
for idx_tmr, tmr in enumerate([0, 5, 10, 15]):
    # Load in histograms and calculate edges
    histograms = np.load('figure5/isi_histograms_tmr_' + str(tmr) + '.npy')
    neural_harm_nums = np.load('figure5/neural_harm_nums_tmr_' + str(tmr) + '.npy')
    edges = np.linspace(start=0, stop=20, num=2200)
    # Plot pcolormesh
    ax[idx_tmr].pcolormesh(edges, neural_harm_nums, histograms, cmap='plasma')
    # Set labels
    ax[idx_tmr].set_xlim((3.5, 4.5))
    # Hide x-axis
    ax[idx_tmr].get_xaxis().set_visible(False)
    # Highlight peaks at F0 multiples
    ax[idx_tmr].plot([4, 4], [4, 12], color=np.array([0, 0, 255])/255, linewidth=2)
    ax[idx_tmr].plot([3/2**(-5.5/12), 3/2**(-5.5/12)], [4, 12], color=np.array([0, 255, 0])/255, linewidth=2)
    ax[idx_tmr].plot([5/2**(6/12), 5/2**(6/12)], [4, 12], color=np.array([0, 255, 255])/255, linewidth=2)
    ax[idx_tmr].plot([6/2**(6/12), 6/2**(6/12)], [4, 12], color=np.array([0, 255, 255])/255, linewidth=2)
# Save to disk
plt.tight_layout()
plt.savefig('plots/fig5a2.png')

# Create subplot (part c)
fig, ax = plt.subplots(4, 1, figsize=(4, 6), sharex='all')
for idx_tmr, tmr in enumerate([0, 5, 10, 15]):
    # Load in histograms and calculate edges
    histograms = np.load('figure5/isi_histograms_tmr_' + str(tmr) + '.npy')
    neural_harm_nums = np.load('figure5/neural_harm_nums_tmr_' + str(tmr) + '.npy')
    edges = np.linspace(start=0, stop=20, num=2200)
    # Plot pcolormesh
    ax[idx_tmr].pcolormesh(edges, neural_harm_nums, histograms, cmap='plasma')
    # Set labels
    ax[idx_tmr].set_xlim((0.5, 1.5))
    # Hide x-axis
    ax[idx_tmr].get_xaxis().set_visible(False)
    # Highlight peaks at F0 multiples
    ax[idx_tmr].plot([1, 1], [4, 12], color=np.array([0, 0, 255])/255, linewidth=2)
    ax[idx_tmr].plot([1/2**(-5.5/12), 1/2**(-5.5/12)], [4, 12], color=np.array([0, 255, 0])/255, linewidth=2)
    ax[idx_tmr].plot([1/2**(6/12), 1/2**(6/12)], [4, 12], color=np.array([0, 255, 255])/255, linewidth=2)
    ax[idx_tmr].plot([2/2**(6/12), 2/2**(6/12)], [4, 12], color=np.array([0, 255, 255])/255, linewidth=2)
# Save to disk
plt.tight_layout()
plt.savefig('plots/fig5a3.png')