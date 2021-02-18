"""
This script estimates firing rates for a pure tone at CF for a range of CFs and then calculates vector strength curves
for those responses.
"""
import apcmodels.synthesis as sy
import apcmodels.simulation as si
import apcmodels.anf as anf
import numpy as np
import config as cfg
import os
from functions import adjust_level


# Define function to estimate vector strength
def estimate_vector_strength(freqs, fs, model, model_name, n_rep):
    """
    Calculates the firing rate of an auditory nerve for a pure tone at a range of freqs (with the CF matched to the
    freq) and then estimates vector strength for those responses. Spike trains are generated from the firing rates
    using a Poisson spike generation multiple times.

    Parameters:
        freqs (ndarray): list or array of frequencies to test, in Hz
        fs (int): sampling rate of the stimulus and model simulation in Hz
        model: model object from apcmodels.anf
        model_name (str): name of model, either 'Heinz2001', 'Zilany2014', or 'Verhulst2018'
        n_rep (int): number of repeated spike trains to generate

    Returns:
        output (ndarray): array of mean firing rates of shape (n_freq, n_level)
    """
    # Write simple wrapper around runfunc to compute Poisson spike trains, calculate vector strength, and average
    def vector_strength_wrapper(runfunc):
        def inner(params):
            # Calculate firing rate
            rates = runfunc(params)  # output is expected to be (1, n_sample)

            # Calculate spike times
            spike_times = [calculate_spike_times(simulate_poisson_spike_train(rates.T, fs), fs) for rep in
                           range(params['n_rep'])]

            # Calculate vector strength for each spike train
            vector_strengths = list(map(lambda x: calculate_vector_strength(x, params['freq']), spike_times))

            # Return mean and standard error
            return np.mean(vector_strengths), np.std(vector_strengths)/np.sqrt(len(vector_strengths))

        def simulate_poisson_spike_train(fr, fs):
            return np.random.rand(len(fr), 1) < fr * (1 / fs)

        def calculate_spike_times(spike_train, fs):
            return np.where(spike_train)[0] * (1 / fs)

        def calculate_vector_strength(spike_train, freq):
            return np.abs(1/len(spike_train) * np.sum(np.exp(1j * 2*np.pi*freq * spike_train)))

        return inner

    # Encode fixed stimulus params in dict
    params = si.Parameters(dur=0.50, dur_ramp=0.020, level=20, phase=0, fs=fs)  # encode fixed params
    params.wiggle_parallel(['freq', 'cf_low', 'cf_high'], [freqs, freqs, freqs])  # wiggle freqs and cfs

    # Add model params to dict
    params.append(['n_cf', 'n_rep'], [1, n_rep])

    # Correct levels to dB re: threshold
    for param in params:
        param['level'] = adjust_level(param['freq'], param['level'], model_name)

    # Synthesize stimuli
    synth = sy.PureTone()
    stimuli = synth.synthesize_sequence(params)
    params.add_inputs(stimuli)

    # Run model
    sim = model()
    if model_name == 'Verhulst2018':
        # If we're using Verhulst model, we should only run one thread as we don't have enough RAM to run default 8
        results = sim.run(params, runfunc=vector_strength_wrapper(sim.simulate), n_thread=1)
    else:
        results = sim.run(params, runfunc=vector_strength_wrapper(sim.simulate))

    # Return average of each response
    return results


# Define parameters
fs = int(200e3)  # sampling rate in Hz
n_rep = 50
freqs = 10**np.linspace(np.log10(200), np.log10(20000), 25)  # CFs for which we will measure vector strength curves

# Save cfs, levels, and freqs to disk
np.save(os.path.join(cfg.root_directory, 'nofigure/vector_strength_curves/freqs.npy'), freqs)

# Loop over models and model names
for model, model_name in zip([anf.AuditoryNerveHeinz2001Numba, anf.AuditoryNerveZilany2014, anf.AuditoryNerveVerhulst2018],
                             ['Heinz2001', 'Zilany2014', 'Verhulst2018']):
    # Compute frequency-level profile for each cf
    vector_strength_curve = estimate_vector_strength(freqs, fs, model, model_name, n_rep)
    # Save neural responses to disk
    np.save(os.path.join(cfg.root_directory, 'nofigure/vector_strength_curves/', model_name + '.npy'),
            vector_strength_curve)



