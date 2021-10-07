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
from util.functions import ISOToneGuest2021_exp1a, ISOToneGuest2021_exp1b, adjust_level


def simulate_supfigure3_f0dls(model, model_name, fs, stimulus, code):
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
    F0s = 10**np.linspace(np.log10(280) - 0.2, np.log10(1400) + 0.1, 24)  # simulate 8th harmonic of F0s
    levels = [20, 30, 40]  # dB SPL, per component
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
    synth = stimulus()
    stimuli = synth.synthesize_sequence(params)
    params.add_inputs(stimuli)

    # Construct simulation and run
    sim = model()
    results = sim.run(params, runfunc=dc.decode_ideal_observer(sim.simulate))

    # Compile results
    save_to_csv([res[0] for res in results], params,
                'supfigure3/' + model_name + '_supfigure3_unroved_AI_' + code + '.csv', decoding_type='AI',
                model=model_name, roving_type='none', code=code)
    save_to_csv([res[1] for res in results], params,
                'supfigure3/' + model_name + '_supfigure3_unroved_RP_' + code + '.csv', decoding_type='RP',
                model=model_name, roving_type='none', code=code)


# Loop through models and calculate FDLs for each model
for stimulus, code in zip([ISOToneGuest2021_exp1a, ISOToneGuest2021_exp1b], ['exp1a', 'exp1b']):
    simulate_supfigure3_f0dls(anf.AuditoryNerveZilany2014, 'Zilany2014', 200e3, stimulus, code)

