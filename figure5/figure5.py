"""
This script implements the simulations described in Figure 5 of Guest and Oxenham (2021).
"""
import apcmodels.simulation as si
import apcmodels.anf as anf
import apcmodels.decode as dc
from apcmodels.util import save_to_csv
import numpy as np
import config as cfg
from functions import ISOToneGuest2021, adjust_level


def simulate_figure5_f0dls(model, model_name, fs):
    """
    Estimates F0 difference limens (FDLs) using ideal observer analysis for a given auditory nerve model. Saves
    the results to disk. This specific harmonic complex tone stimulus used in this simulation is from Guest and
    Oxenham (2021), although no acoustic noise is included in the stimulus.

    Arguments:
        model: model object from apcmodels.anf

        model_name (str): name of the model

        fs (int): sampling rate in Hz
    """
    # Define stimulus parameters
    F0s = 10**np.linspace(np.log10(280) - 0.2, np.log10(1400) + 0.1, 24)  # simulate 8th harmonic of F0s
    levels = [20, 30, 40]  # dB SPL, per component
    dur = 0.10  # seconds
    dur_ramp = 0.01  # seconds

    # Define model parameters
    cf_low = 5*F0s
    cf_high = 11*F0s
    n_cf = 40
    n_fiber_per_chan = 40  # TODO: check this value

    # Encode parameters
    params = {'dur': dur, 'dur_ramp': dur_ramp, 'fs': fs}  # encode fixed parameters in a dict
    params = si.wiggle_parameters(params, 'F0', F0s)  # wiggle frequencies
    params = si.stitch_parameters(params, 'cf_low', cf_low)  # stitch cf_low (each F0 corresponds to a cf_low)
    params = si.stitch_parameters(params, 'cf_high', cf_high)  # stitch cf_high (each F0 corresponds to a cf_high)
    params = si.wiggle_parameters(params, 'level', levels)  # wiggle levels
    params = si.append_parameters(params, ['n_cf', 'delta_theta', 'API', 'n_fiber_per_chan'],
                                  [n_cf, [0.001], np.zeros(1), n_fiber_per_chan])  # append other model parameters

    # Adjust levels
    params = si.flatten_parameters(params)  # flatten out params
    for ele in params:
        ele['nominal_level'] = ele['level']
        ele['level'] = adjust_level(ele['F0']*8, ele['level'], model_name)

    # Encode increments
    params = si.increment_parameters(params, {'F0': 0.001})  # increment frequency

    # Synthesize stimuli
    synth = ISOToneGuest2021()
    stimuli = synth.synthesize_sequence(params)
    params = si.stitch_parameters(params, '_input', stimuli)

    # Construct simulation and run
    sim = model()
    results = sim.run(params, runfunc=dc.decode_ideal_observer(sim.simulate), parallel=True, progress=True)

    # Compile results
    save_to_csv([res[0] for res in results], params,
                cfg.root_directory + 'figure5/' + model_name + '_figure5_unroved_AI.csv', decoding_type='AI',
                model=model_name, roving_type='none')
    save_to_csv([res[1] for res in results], params,
                cfg.root_directory + 'figure5/' + model_name + '_figure5_unroved_RP.csv', decoding_type='RP',
                model=model_name, roving_type='none')


# Loop through models and calculate FDLs for each model
for model, model_name, fs in zip([anf.AuditoryNerveHeinz2001Numba],
                                 ['Heinz2001'],
                                 [int(1000e3)]):
    simulate_figure5_f0dls(model, model_name, fs)
