import apcmodels.simulation as si
import apcmodels.anf as anf
import numpy as np
import itertools
import os, sys
sys.path.append(os.getcwd())
from util.functions import ComplexToneLarsen2008


def calculate_ISI_histograms_Larsen_2008(cf, neural_harm_nums):
    """
    Synthesizes a complex tone dyad and simulates interspike intervals for an auditory nerve fiber for it. This is
    performed at a range of neural harmonic numbers.

    Parameters:
        cf (float): cf of the complex tone
        neural_harm_nums (np.ndarray): array of neural harmonic numbers (CF/F0) to test

    Returns:
        histograms (list): list of histograms counts for ISI histogram for each neural harmonic number
    """
    # Setup params
    params = si.Parameters(fs=int(200e3), n_cf=1, cf_low=cf, cf_high=cf, level=70, anf_num=(0, 1, 0), species='cat')
    params.wiggle('F0', cf / neural_harm_nums)

    # Synthesize stimuli and add to params
    stim = ComplexToneLarsen2008()
    params.add_inputs(stim.synthesize_sequence(params))

    # Encode repeats
    params.repeat(200)

    # Select model and run
    sim = anf.AuditoryNerveZilany2014()
    results = sim.run(params, runfunc=lambda x: [sim.simulate_spikes(ele) for ele in x])

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
        ISIs_period = np.array(ISIs) / (1 / (cf / neural_harm_nums[idx_result]))
        # Calculate histogram and append to results
        hist_low, edges_low = np.histogram(ISIs_period, bins=2200, range=(0, 20))
        histograms.append(hist_low)

    # Return
    return histograms


# Calculate ISI histograms for one Larsen and Delgutte (2008) stimulus and save to disk
histograms = calculate_ISI_histograms_Larsen_2008(816, np.linspace(2, 6, 40))
np.save('supplemental_figure_1/figure7e/isi_histograms.npy', histograms)
np.save('supplemental_figure_1/figure7e/neural_harm_nums.npy', np.linspace(2, 6, 40))

