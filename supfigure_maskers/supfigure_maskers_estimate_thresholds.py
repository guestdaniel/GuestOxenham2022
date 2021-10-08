"""
This script implements the simulations described in Figure 6 of Guest and Oxenham (2021).
"""
from re import finditer
import apcmodels.simulation as si
import apcmodels.anf as anf
import apcmodels.decode as dc
from apcmodels.util import save_to_csv
import numpy as np
import os, sys
sys.path.append(os.getcwd())
from util.functions import ISOToneGuest2021, GEOMToneGuest2021, adjust_level


def simulate(model, model_name, fs, stim='iso', delta=0.001, finite_difference_method='forward', masker_interval=1):
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
    if finite_difference_method == 'forward' or finite_difference_method == 'backward':
        F0s = np.array([280, 1400])
        nominal_F0s = np.array([280, 1400])
    else:
        F0s = np.array([280, 1400]) - delta
        nominal_F0s = np.array([280, 1400])
    F0s_masker = F0s*2**(masker_interval/12)
    levels = [30]  # dB SPL, per component
    dur = 0.10  # seconds
    dur_ramp = 0.01  # seconds

    # Define model parameters
    cf_low = 5*F0s
    cf_high = 11*F0s
    n_cf = 40
    n_fiber_per_chan = round(((np.log10(11/5)/3)*18000)/n_cf)  # assume ~18k HSR fibers from 0.2 to 20 kHz

    # Encode parameters
    params = si.Parameters(dur=dur, dur_ramp=dur_ramp, fs=fs, n_cf=n_cf, API=np.zeros(1),
                           n_fiber_per_chan=n_fiber_per_chan, model_name=model_name)
    params.wiggle('F0', F0s)                                   # wiggle F0s
    params.stitch('nominal_F0', nominal_F0s)
    params.stitch('F0_masker', F0s_masker)
    params.stitch('cf_low', cf_low)                            # stitch cf_low (each freq corresponds to a cf_low)
    params.stitch('cf_high', cf_high)                          # stitch cf_high (each freq corresponds to a cf_high)
    params.wiggle('level', levels)                             # wiggle levels
    if finite_difference_method == 'forward':
        params.append('delta_theta', [delta])
    elif finite_difference_method == 'backward':
        params.append('delta_theta', [delta])
    elif finite_difference_method == 'central':
        params.append('delta_theta', [delta*2])

    # Adjust levels to be in dB re: threshold
    params.flatten()
    for ele in params:
        ele['nominal_level'] = ele['level']                                 # encode nominal level (dB re: threshold)
        ele['level'] = adjust_level(ele['F0']*8, ele['level'], model_name)  # encode actual level (dB SPL)

    # Encode increments
    if finite_difference_method == 'forward':
        params.increment({'F0': delta})  # increment F0
    elif finite_difference_method == 'backward':
        params.increment({'F0': -1*delta})
    elif finite_difference_method == 'central':
        params.increment({'F0': 2*delta})

    # Synthesize stimuli
    if stim == 'iso':
        synth = ISOToneGuest2021()
    else:
        synth = GEOMToneGuest2021()
    stimuli = synth.synthesize_sequence(params)
    params.add_inputs(stimuli)

    # Construct simulation and run
    sim = model()
    results = sim.run(params, runfunc=dc.decode_ideal_observer(sim.simulate))

    # Compile results
    save_to_csv([res[0] for res in results], params,
                'supfigure_maskers/' + model_name + '_supfigure_maskers_unroved_AI_' + str(delta) + '_' + finite_difference_method + '_' + stim + '.csv', decoding_type='AI',
                model=model_name, roving_type='none', stimulus=stim, delta=delta, finite_difference_method=finite_difference_method)
    save_to_csv([res[1] for res in results], params,
                'supfigure_maskers/' + model_name + '_supfigure_maskers_unroved_RP_' + str(delta) + '_' + finite_difference_method + '_' + stim + '.csv', decoding_type='RP',
                model=model_name, roving_type='none', stimulus=stim, delta=delta, finite_difference_method=finite_difference_method)


# Loop through models and calculate FDLs for each model
for model, model_name, fs in zip([anf.AuditoryNerveHeinz2001],
                                 ['Heinz2001'],
                                 [int(1000e3)]):
    for delta in [1e-6, 1e-4, 1e-2, 1e0, 1e1]:
        for finite_differences_method in ['forward', 'backward', 'central']:
            for stim in ['iso', 'geom']:
                simulate(model, model_name, fs, delta=delta, stim=stim, 
                         finite_difference_method=finite_differences_method, masker_interval=0.5)
