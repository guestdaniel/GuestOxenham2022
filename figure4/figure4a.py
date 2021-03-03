"""
This script implements the simulations described in Figure 4 of Guest and Oxenham (2021).
"""
import apcmodels.synthesis as sy
import apcmodels.simulation as si
import apcmodels.anf as anf
import apcmodels.decode as dc
from apcmodels.util import save_to_csv
import numpy as np
import os, sys
sys.path.append(os.getcwd())
from util.functions import adjust_level


def simulate_figure4_fdls(model, model_name, fs):
    """
    Estimates frequency difference limens (FDLs) using ideal observer analysis for a given auditory nerve model. Saves
    the results to disk.

    Arguments:
        model: model object from apcmodels.anf

        model_name (str): name of the model

        fs (int): sampling rate in Hz
    """
    # Define stimulus parameters
    freqs = 8 * 10**np.linspace(np.log10(280) - 0.2, np.log10(1400) + 0.1, 24)  # simulate 8th harmonic of F0s
    levels = [20, 30, 40]  # dB SPL
    dur = 0.10  # seconds
    dur_ramp = 0.01  # seconds

    # Define model parameters
    cf_low = 0.5*freqs
    cf_high = 1.5*freqs
    n_cf = 40
    n_fiber_per_chan = 40  # TODO: check this value

    # Encode parameters
    params = {'dur': dur, 'dur_ramp': dur_ramp, 'fs': fs}  # encode fixed parameters in a dict
    params = si.wiggle_parameters(params, 'freq', freqs)  # wiggle frequencies
    params = si.stitch_parameters(params, 'cf_low', cf_low)  # stitch cf_low (each frequency corresponds to a cf_low)
    params = si.stitch_parameters(params, 'cf_high', cf_high)  # stitch cf_high (each frequency corresponds to a cf_high)
    params = si.wiggle_parameters(params, 'level', levels)  # wiggle levels
    params = si.append_parameters(params, ['n_cf', 'delta_theta', 'API', 'n_fiber_per_chan', 'model_name'],
                                  [n_cf, [0.001], np.zeros(1), n_fiber_per_chan, model_name])  # append other model parameters

    # Adjust levels to be in dB re: threshold
    params = si.flatten_parameters(params)  # flatten out params
    for ele in params:
        ele['nominal_level'] = ele['level']  # save the nominal level (dB re: threshold)
        ele['level'] = adjust_level(ele['freq'], ele['level'], model_name)  # encode the actual level (dB SPL)

    # Encode increments
    params = si.increment_parameters(params, {'freq': 0.001})  # increment frequency

    # Synthesize stimuli
    synth = sy.PureTone()
    stimuli = synth.synthesize_sequence(params)
    params = si.stitch_parameters(params, '_input', stimuli)

    # Construct simulation and run
    sim = model()
    results = sim.run(params, runfunc=dc.decode_ideal_observer(sim.simulate), parallel=True, hide_progress=False)

    # Compile results
    save_to_csv([res[0] for res in results], params,
                'figure4/' + model_name + '_figure4_unroved_AI.csv', decoding_type='AI',
                model=model_name, roving_type='none')
    save_to_csv([res[1] for res in results], params,
                'figure4/' + model_name + '_figure4_unroved_RP.csv', decoding_type='RP',
                model=model_name, roving_type='none')


# Loop through models and calculate FDLs for each model
for model, model_name, fs in zip([anf.AuditoryNerveHeinz2001Numba, anf.AuditoryNerveZilany2014, anf.AuditoryNerveVerhulst2018],
                                 ['Heinz2001', 'Zilany2014', 'Verhulst2018'],
                                 [int(1000e3), int(200e3), int(200e3)]):
    simulate_figure4_fdls(model, model_name, fs)
