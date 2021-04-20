"""
This script implements the simulations described in Figure 6b of Guest and Oxenham (2021).
"""
import apcmodels.simulation as si
import apcmodels.anf as anf
import apcmodels.decode as dc
from apcmodels.util import save_to_csv
import numpy as np
import os, sys
sys.path.append(os.getcwd())
from util.functions import adjust_level, ISOToneGuest2021_exp1a


def simulate_figure7_f0dls_phase_roving(model, model_name, fs, n_rep=10):
    """
    Estimates F0 difference limens (F0DLs) using ideal observer analysis for a given auditory nerve model. Saves
    the results to disk. The simulations include phase randomization (specifically, each component from 6-10 has its
    starting phase randomized on each trial).

    Args:
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
    n_fiber_per_chan = round(((np.log10(11/5)/3)*18000)/n_cf)  # assume ~18k HSR fibers from 0.2 to 20 kHz

    # Construct delta_theta and API matrices/arrays
    delta_theta = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001]
    API = np.zeros((6, 6))
    API[np.diag_indices(6)] = [0, 1/360**2, 1/360**2, 1/360**2, 1/360**2, 1/360**2]

    # Encode parameters
    params = si.Parameters(dur=dur, dur_ramp=dur_ramp, fs=fs, n_cf=n_cf, delta_theta=delta_theta, API=API,
                           n_fiber_per_chan=n_fiber_per_chan, model_name=model_name)
    params.append('phase', lambda: np.random.uniform(0, 360, 5))  # append random phase
    params.wiggle('F0', F0s)                                   # wiggle F0s
    params.stitch('cf_low', cf_low)                            # stitch cf_low (each freq corresponds to a cf_low)
    params.stitch('cf_high', cf_high)                          # stitch cf_high (each freq corresponds to a cf_high)
    params.wiggle('level', levels)                             # wiggle levels

    # Adjust levels to be in dB re: threshold
    params.flatten()
    for ele in params:
        ele['nominal_level'] = ele['level']                                 # encode nominal level (dB re: threshold)
        ele['level'] = adjust_level(ele['F0']*8, ele['level'], model_name)  # encode actual level (dB SPL)

    # Encode repeats and increments
    params.repeat(n_rep)
    params.increment({'F0': 0.001,                                  # increment F0
                      '(1)_phase': np.array([0.001, 0, 0, 0, 0]),   # increment phase of H6
                      '(2)_phase': np.array([0, 0.001, 0, 0, 0]),   # increment phase of H7
                      '(3)_phase': np.array([0, 0, 0.001, 0, 0]),   # increment phase of H8
                      '(4)_phase': np.array([0, 0, 0, 0.001, 0]),   # increment phase of H9
                      '(5)_phase': np.array([0, 0, 0, 0, 0.001])})  # increment phase of H10

    # Synthesize stimuli
    synth = ISOToneGuest2021_exp1a()
    stimuli = synth.synthesize_sequence(params)
    params.add_inputs(stimuli)

    # Construct simulation and run
    sim = model()
    results = sim.run(params, runfunc=dc.decode_ideal_observer(sim.simulate))

    # Compile results
    save_to_csv([res[0] for res in results], params,
                'figure7/' + model_name + '_figure7_phase_roving_AI.csv', decoding_type='AI',
                model=model_name, roving_type='phase')
    save_to_csv([res[1] for res in results], params,
                'figure7/' + model_name + '_figure7_phase_roving_RP.csv', decoding_type='RP',
                model=model_name, roving_type='phase')


def simulate_figure7_f0dls_level_roving(model, model_name, fs, n_rep=10):
    """
    Estimates F0 difference limens (F0DLs) using ideal observer analysis for a given auditory nerve model. Saves
    the results to disk. Includes simulations of level roving.

    Args:
        model: model object from apcmodels.anf
        model_name (str): name of the model
        fs (int): sampling rate in Hz
        n_rep (int): number of repetitions to perform for each simulation
    """
    # Define stimulus parameters
    F0s = 10**np.linspace(np.log10(280) - 0.2, np.log10(1400) + 0.1, 24)  # simulate 8th harmonic of F0s
    levels = [lambda: np.random.uniform(17, 23, 5),
              lambda: np.random.uniform(27, 33, 5),
              lambda: np.random.uniform(37, 43, 5)]  # dB SL
    nominal_levels = [20, 30, 40]  # dB SL
    dur = 0.10  # seconds
    dur_ramp = 0.01  # seconds

    # Define model parameters
    cf_low = 5*F0s
    cf_high = 11*F0s
    n_cf = 40
    n_fiber_per_chan = round(((np.log10(11/5)/3)*18000)/n_cf)  # assume ~18k HSR fibers from 0.2 to 20 kHz

    # Construct delta_theta and API matrices/arrays
    delta_theta = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001]
    API = np.zeros((6, 6))
    API[np.diag_indices(6)] = [0, 1/6**2, 1/6**2, 1/6**2, 1/6**2, 1/6**2]

    # Encode parameters
    params = si.Parameters(dur=dur, dur_ramp=dur_ramp, fs=fs, n_cf=n_cf, delta_theta=delta_theta, API=API,
                           n_fiber_per_chan=n_fiber_per_chan, model_name=model_name)
    params.wiggle('F0', F0s)                                   # wiggle F0s
    params.stitch('cf_low', cf_low)                            # stitch cf_low (each freq corresponds to a cf_low)
    params.stitch('cf_high', cf_high)                          # stitch cf_high (each freq corresponds to a cf_high)
    params.wiggle_parallel(['level', 'nominal_level'], [levels, nominal_levels])  # wiggle levels

    # Flatten parameters
    params.flatten()

    # Encode repeats and increments
    params.repeat(n_rep)
    params.increment({'F0': 0.001,                                  # increment F0
                      '(1)_level': np.array([0.001, 0, 0, 0, 0]),   # increment level of H6
                      '(2)_level': np.array([0, 0.001, 0, 0, 0]),   # increment level of H7
                      '(3)_level': np.array([0, 0, 0.001, 0, 0]),   # increment level of H8
                      '(4)_level': np.array([0, 0, 0, 0.001, 0]),   # increment level of H9
                      '(5)_level': np.array([0, 0, 0, 0, 0.001])})  # increment level of H10

    # Adjust levels to be in dB re: threshold
    for ele in params:   # loop through elements of params
        for repeat in ele:    # loop through repeats of params
            for inc in repeat:  # loop through increments of params
                inc['level'] = adjust_level(inc['F0']*8, inc['level'], model_name)  # encode the actual level (dB SPL)

    # Synthesize stimuli
    synth = ISOToneGuest2021_exp1a()
    stimuli = synth.synthesize_sequence(params)
    params.add_inputs(stimuli)

    # Construct simulation and run
    sim = model()
    results = sim.run(params, runfunc=dc.decode_ideal_observer(sim.simulate))

    # Compile results
    save_to_csv([res[0] for res in results], params,
                'figure7/' + model_name + '_figure7_level_roving_AI.csv', decoding_type='AI',
                model=model_name, roving_type='level')
    save_to_csv([res[1] for res in results], params,
                'figure7/' + model_name + '_figure7_level_roving_RP.csv', decoding_type='RP',
                model=model_name, roving_type='level')


# Loop through models and calculate FDLs for each model (for this figure, we only use Zilany et al. [2014])
for model, model_name, fs in zip([anf.AuditoryNerveZilany2014],
                                 ['Zilany2014'],
                                 [int(200e3)]):
    simulate_figure7_f0dls_phase_roving(model, model_name, fs)
    simulate_figure7_f0dls_level_roving(model, model_name, fs)

