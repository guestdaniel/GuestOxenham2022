"""
This script implements the simulations described in Figure 6 of Guest and Oxenham (2021).
"""
import apcmodels.simulation as si
import apcmodels.anf as anf
import apcmodels.decode as dc
from apcmodels.util import save_to_csv
import numpy as np
import os, sys
sys.path.append(os.getcwd())
from util.functions import ISOToneGuest2021, GEOMToneGuest2021, adjust_level


def simulate_iso(model, model_name, fs):
    """
    Estimates F0 difference limens (FDLs) using ideal observer analysis for a given auditory nerve model. Saves
    the results to disk. This specific harmonic complex tone stimulus used in this simulation is from Guest and
    Oxenham (2021), although no acoustic noise is included in the stimulus.

    Args:
        model: model object from apcmodels.anf
        model_name (str): name of the model
        fs (int): sampling rate in Hz
    """
    # Define stimulus parameters
    #F0s = 10**np.linspace(np.log10(280) - 0.2, np.log10(1400) + 0.1, 24)  # simulate 8th harmonic of F0s
    F0s = np.array([280, 1400])
    levels = [10, 20, 30, 40]  # dB SPL, per component
    dur = 0.10  # seconds
    dur_ramp = 0.01  # seconds

    # Define model parameters
    cf_low = 5*F0s
    cf_high = 11*F0s
    n_cf = 40
    n_fiber_per_chan = round(((np.log10(11/5)/3)*18000)/n_cf)  # assume ~18k HSR fibers from 0.2 to 20 kHz

    # Encode parameters
    params = si.Parameters(dur=dur, dur_ramp=dur_ramp, fs=fs, n_cf=n_cf, delta_theta=[0.001], API=np.zeros(1),
                           n_fiber_per_chan=n_fiber_per_chan, model_name=model_name)
    params.wiggle('F0', F0s)                                   # wiggle F0s
    params.stitch('nominal_F0', F0s)
    params.stitch('cf_low', cf_low)                            # stitch cf_low (each freq corresponds to a cf_low)
    params.stitch('cf_high', cf_high)                          # stitch cf_high (each freq corresponds to a cf_high)
    params.wiggle('level', levels)                             # wiggle levels

    # Adjust levels to be in dB re: threshold
    params.flatten()
    for ele in params:
        ele['nominal_level'] = ele['level']                                 # encode nominal level (dB re: threshold)
        ele['level'] = adjust_level(ele['F0']*8, ele['level'], model_name)  # encode actual level (dB SPL)

    # Encode increments
    params.increment({'F0': 0.001})  # increment F0

    # Synthesize stimuli
    synth = ISOToneGuest2021()
    stimuli = synth.synthesize_sequence(params)
    params.add_inputs(stimuli)

    # Construct simulation and run
    sim = model()
    results = sim.run(params, runfunc=dc.decode_ideal_observer(sim.simulate))

    # Compile results
    save_to_csv([res[0] for res in results], params,
                'supfigure_maskers/' + model_name + '_supfigure_maskers_unroved_AI_iso.csv', decoding_type='AI',
                model=model_name, roving_type='none', stimulus='iso')
    save_to_csv([res[1] for res in results], params,
                'supfigure_maskers/' + model_name + '_supfigure_maskers_unroved_RP_iso.csv', decoding_type='RP',
                model=model_name, roving_type='none', stimulus='iso')


def simulate_geom(model, model_name, fs, n_rep=10):
    """
    Estimates F0 difference limens (FDLs) using ideal observer analysis for a given auditory nerve model. Saves
    the results to disk. This specific harmonic complex tone stimulus used in this simulation is from Guest and
    Oxenham (2021), although no acoustic noise is included in the stimulus.

    Args:
        model: model object from apcmodels.anf
        model_name (str): name of the model
        fs (int): sampling rate in Hz
    """
    # Define stimulus parameters
    #F0s = 10**np.linspace(np.log10(280) - 0.2, np.log10(1400) + 0.1, 24)  # simulate 8th harmonic of F0s
    F0s = np.array([lambda: 280*np.random.uniform(0.9, 1.1, 1)[0], lambda: 1400*np.random.uniform(0.9, 1.1, 1)[0]])
    F0s_masker = np.array([lambda: 280*np.random.uniform(0.9, 1.1, 1)[0], lambda: 1400*np.random.uniform(0.9, 1.1, 1)[0]])
    levels = [10, 20, 30, 40]  # dB SPL, per component
    dur = 0.10  # seconds
    dur_ramp = 0.01  # seconds

    # Define model parameters
    cf_low = np.array([5*280, 5*1400])
    cf_high = np.array([11*280, 11*1400])
    n_cf = 40
    n_fiber_per_chan = round(((np.log10(11/5)/3)*18000)/n_cf)  # assume ~18k HSR fibers from 0.2 to 20 kHz

    # Encode parameters
    params = si.Parameters(dur=dur, dur_ramp=dur_ramp, fs=fs, n_cf=n_cf, delta_theta=[0.001, 0.001], API=np.array([[0, 0], [0, 0]]),
                           n_fiber_per_chan=n_fiber_per_chan, model_name=model_name)
    params.wiggle('F0', F0s)                                   # wiggle F0s
    params.stitch('nominal_F0', np.array([280, 1400]))
    params.stitch('F0_masker', F0s_masker)
    params.stitch('cf_low', cf_low)                            # stitch cf_low (each freq corresponds to a cf_low)
    params.stitch('cf_high', cf_high)                          # stitch cf_high (each freq corresponds to a cf_high)
    params.wiggle('level', levels)                             # wiggle levels

    # Adjust levels to be in dB re: threshold
    params.flatten()
    for ele in params:
        ele['nominal_level'] = ele['level']                                 # encode nominal level (dB re: threshold)
        ele['level'] = adjust_level(ele['nominal_F0']*8, ele['level'], model_name)  # encode actual level (dB SPL)

    # Encode increments
    params.repeat(n_rep)
    params.increment({'F0': 0.001, 'F0_masker': 0.001})  # increment F0

    # Synthesize stimuli
    synth = GEOMToneGuest2021()
    stimuli = synth.synthesize_sequence(params)
    params.add_inputs(stimuli)

    # Construct simulation and run
    sim = model()
    results = sim.run(params, runfunc=dc.decode_ideal_observer(sim.simulate))

    # Compile results
    save_to_csv([res[0] for res in results], params,
                'supfigure_maskers/' + model_name + '_supfigure_maskers_unroved_AI_geom.csv', decoding_type='AI',
                model=model_name, roving_type='none', stimulus='geom')
    save_to_csv([res[1] for res in results], params,
                'supfigure_maskers/' + model_name + '_supfigure_maskers_unroved_RP_geom.csv', decoding_type='RP',
                model=model_name, roving_type='none', stimulus='geom')


# Loop through models and calculate FDLs for each model
for model, model_name, fs in zip([anf.AuditoryNerveHeinz2001, anf.AuditoryNerveZilany2014],
                                 ['Heinz2001', 'Zilany2014'],
                                 [int(1250e3), int(500e3)]):
    simulate_iso(model, model_name, fs)
    simulate_geom(model, model_name, fs)

