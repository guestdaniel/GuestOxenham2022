"""
This script implements the simulations described in Figure 4b of Guest and Oxenham (2021).
"""
import apcmodels.simulation as si
import apcmodels.anf as anf
import apcmodels.decode as dc
from apcmodels.util import save_to_csv
import numpy as np
import config as cfg
from functions import adjust_level, ISOToneGuest2021_exp1a


def simulate_figure5_f0dls_phase_roving(model, model_name, fs, n_rep=10):
    """
    Estimates F0 difference limens (F0DLs) using ideal observer analysis for a given auditory nerve model. Saves
    the results to disk. The simulations include phase randomization (specifically, each component from 6-10 has its
    starting phase randomized on each trial).

    Arguments:
        model: model object from apcmodels.anf

        model_name (str): name of the model

        fs (int): sampling rate in Hz

        n_rep (int): number of repetitions to perform for each simulation
    """
    # Define stimulus parameters
    F0s = 10**np.linspace(np.log10(280) - 0.2, np.log10(1400) + 0.1, 24)  # simulate 8th harmonic of F0s
    levels = [20, 30, 40]  # dB SPL
    dur = 0.10  # seconds
    dur_ramp = 0.01  # seconds

    # Define model parameters
    cf_low = 5*F0s
    cf_high = 11*F0s
    n_cf = 40
    n_fiber_per_chan = 40  # TODO: check this value

    # Encode parameters
    delta_theta = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001]
    API = np.zeros((6, 6))
    API[np.diag_indices(6)] = [0, 1/360**2, 1/360**2, 1/360**2, 1/360**2, 1/360**2]
    params = {'dur': dur, 'dur_ramp': dur_ramp, 'fs': fs}  # encode fixed parameters in a dict
    params['phase'] = lambda: np.random.uniform(0, 360, 5)
    params = si.wiggle_parameters(params, 'F0', F0s)  # wiggle frequencies
    params = si.stitch_parameters(params, 'cf_low', cf_low)  # stitch cf_low (each F0 corresponds to a cf_low)
    params = si.stitch_parameters(params, 'cf_high', cf_high)  # stitch cf_high (each F0 corresponds to a cf_high)
    params = si.wiggle_parameters(params, 'level', levels)  # wiggle levels
    params = si.append_parameters(params, ['n_cf', 'delta_theta', 'API', 'n_fiber_per_chan', 'model_name'],
                                  [n_cf, delta_theta, API, n_fiber_per_chan, model_name])  # append other model parameters

    # Adjust levels to be in dB re: threshold
    params = si.flatten_parameters(params)  # flatten out params
    for ele in params:
        ele['nominal_level'] = ele['level']  # save the nominal level (dB re: threshold)
        ele['level'] = adjust_level(ele['F0']*8, ele['level'], model_name)  # encode the actual level (dB SPL)

    # Encode repeats and increments
    params = si.repeat_parameters(params, n_rep)
    params = si.increment_parameters(params, {'F0': 0.001,                                  # increment F0
                                              '(1)_phase': np.array([0.001, 0, 0, 0, 0]),   # increment phase of H6
                                              '(2)_phase': np.array([0, 0.001, 0, 0, 0]),   # increment phase of H7
                                              '(3)_phase': np.array([0, 0, 0.001, 0, 0]),   # increment phase of H8
                                              '(4)_phase': np.array([0, 0, 0, 0.001, 0]),   # increment phase of H9
                                              '(5)_phase': np.array([0, 0, 0, 0, 0.001])})  # increment phase of H10

    # Synthesize stimuli
    synth = ISOToneGuest2021_exp1a()
    stimuli = synth.synthesize_sequence(params)
    params = si.stitch_parameters(params, '_input', stimuli)

    # Construct simulation and run
    sim = model()
    results = sim.run(params, runfunc=dc.decode_ideal_observer(sim.simulate), parallel=True, hide_progress=False)

    # Compile results
    save_to_csv([res[0] for res in results], params,
                'figure5/' + model_name + '_figure5_phase_roving_AI.csv', decoding_type='AI',
                model=model_name, roving_type='phase')
    save_to_csv([res[1] for res in results], params,
                'figure5/' + model_name + '_figure5_phase_roving_RP.csv', decoding_type='RP',
                model=model_name, roving_type='phase')


def simulate_figure5_f0dls_level_roving(model, model_name, fs, n_rep=10):
    """
    Estimates F0 difference limens (F0DLs) using ideal observer analysis for a given auditory nerve model. Saves
    the results to disk. Includes simulations of level roving.

    Arguments:
        model: model object from apcmodels.anf

        model_name (str): name of the model

        fs (int): sampling rate in Hz

        n_rep (int): number of repetitions to perform for each simulation
    """
    # Define stimulus parameters
    F0s = 10**np.linspace(np.log10(280) - 0.2, np.log10(1400) + 0.1, 24)  # simulate 8th harmonic of F0s
    levels = [np.random.uniform(17, 23, 5), np.random.uniform(27, 33, 5), np.random.uniform(37, 43, 5)]  # dB SL
    nominal_levels = [20, 30, 40]
    dur = 0.10  # seconds
    dur_ramp = 0.01  # seconds

    # Define model parameters
    cf_low = 5*F0s
    cf_high = 11*F0s
    n_cf = 40
    n_fiber_per_chan = 40  # TODO: check this value

    # Encode parameters
    delta_theta = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001]
    API = np.zeros((6, 6))
    API[np.diag_indices(6)] = [0, 1/6**2, 1/6**2, 1/6**2, 1/6**2, 1/6**2]
    params = {'dur': dur, 'dur_ramp': dur_ramp, 'fs': fs}  # encode fixed parameters in a dict
    params = si.wiggle_parameters(params, 'F0', F0s)  # wiggle frequencies
    params = si.stitch_parameters(params, 'cf_low', cf_low)  # stitch cf_low (each F0 corresponds to a cf_low)
    params = si.stitch_parameters(params, 'cf_high', cf_high)  # stitch cf_high (each F0 corresponds to a cf_high)
    params = si.wiggle_parameters_parallel(params, ['level', 'nominal_level'], [levels, nominal_levels])  # wiggle levels
    params = si.append_parameters(params, ['n_cf', 'delta_theta', 'API', 'n_fiber_per_chan', 'model_name'],
                                  [n_cf, delta_theta, API, n_fiber_per_chan, model_name])  # append other model parameters

    # Adjust levels to be in dB re: threshold
    params = si.flatten_parameters(params)  # flatten out params

    # Encode repeats and increments
    params = si.repeat_parameters(params, n_rep)
    params = si.increment_parameters(params, {'F0': 0.001,                                  # increment F0
                                              '(1)_level': np.array([0.001, 0, 0, 0, 0]),   # increment level of H6
                                              '(2)_level': np.array([0, 0.001, 0, 0, 0]),   # increment level of H7
                                              '(3)_level': np.array([0, 0, 0.001, 0, 0]),   # increment level of H8
                                              '(4)_level': np.array([0, 0, 0, 0.001, 0]),   # increment level of H9
                                              '(5)_level': np.array([0, 0, 0, 0, 0.001])})  # increment level of H10

    # Adjust levels to be in dB re: threshold
    for ele in params:   # loop through elements of params
        for x in ele:    # loop through repeats of params
            for y in x:  # loop through increments of params
                y['level'] = adjust_level(y['F0']*8, y['level'], model_name)  # encode the actual level (dB SPL)

    # Synthesize stimuli
    synth = ISOToneGuest2021_exp1a()
    stimuli = synth.synthesize_sequence(params)
    params = si.stitch_parameters(params, '_input', stimuli)

    # Construct simulation and run
    sim = model()
    results = sim.run(params, runfunc=dc.decode_ideal_observer(sim.simulate), parallel=True, hide_progress=False)

    # Compile results
    save_to_csv([res[0] for res in results], params,
                'figure5/' + model_name + '_figure5_level_roving_AI.csv', decoding_type='AI',
                model=model_name, roving_type='level')
    save_to_csv([res[1] for res in results], params,
                'figure5/' + model_name + '_figure5_level_roving_RP.csv', decoding_type='RP',
                model=model_name, roving_type='level')


# Loop through models and calculate FDLs for each model
for model, model_name, fs in zip([anf.AuditoryNerveHeinz2001Numba],
                                 ['Heinz2001'],  # TODO: switch to Zilany
                                 [int(500e3)]):
    simulate_figure5_f0dls_phase_roving(model, model_name, fs)
    simulate_figure5_f0dls_level_roving(model, model_name, fs)

