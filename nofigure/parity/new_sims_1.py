"""
This script implements the simulations described in Figure 6 of Guest and Oxenham (2021).
"""
import apcmodels.simulation as si
import apcmodels.anf as anf
import apcmodels.decode as dc
import apcmodels.synthesis as sy
import apcmodels.signal as sg
from apcmodels.util import save_to_csv
import numpy as np
import os, sys
sys.path.append(os.getcwd())
from util.functions import ISOToneGuest2021_exp1a, adjust_level

class ISOToneGuest2021_exp1a(sy.Synthesizer):
    """ Synthesizes the ISO stimulus in Guest and Oxenham (2021).

    Simplified version of the stimulus in Guest and Oxenham (2021). This version is the version from Experiment 1a, and
    does not include acoustic masking noise.
    """
    def __init__(self):
        super().__init__(stimulus_name='ISO Tone')

    def synthesize(self, F0, comps, dur=0.350, dur_ramp=0.02, level=None, phase=None, fs=int(48e3), **kwargs):
        """
        Synthesizes a harmonic complex tone composed of components 6-10 of the F0. This is the same stimulus as used
        in Experiment 1a. An additional argument, comps, allows the user to specify which harmonics to synthesize. 

        Arguments:
            F0 (float): F0 of the complex tone in Hz
            level (float, ndarray): level per-component of the complex tone in dB SPL, can be either a float (in this
                case, the same level is used for all components) or an ndarray indicating the level of each component
            phase (float, ndarray): phase offset applied to each component of the complex tone in degrees, can be
                either a float (in which case the same phase offset is used for all components) or an ndarray indicating
                the phase offset of each component.
            dur (float): duration in seconds
            dur_ramp (float): duration of raised-cosine ramp in seconds
            fs (int): sampling rate in Hz

        Returns:
            output (array): complex tone stimulus
        """
        # Create array of frequencies, levels, and phases
        freqs = F0*comps
        if level is None:
            level = 40*np.ones(len(freqs))  # default to 40 dB SPL per component
        elif isinstance(level, float) or isinstance(level, int) is int:
            level = level*np.ones(len(freqs))  # default to 40 dB SPL per component
        if phase is None:
            phase = np.zeros(len(freqs))  # default to sine phase
        elif isinstance(phase, float) or isinstance(phase, int) is int:
            phase = phase + np.zeros(len(freqs))
        # Synthesize, filter, and ramp complex tone signal
        signal = sg.complex_tone(freqs, level, phase, dur, fs)
        signal = sg.cosine_ramp(signal, dur_ramp, fs)
        # Return
        return signal


def simulate_f0dls_dichotic(model, model_name, fs, comps, tag):
    """
    Estimates F0 difference limens (FDLs) using ideal observer analysis for a given auditory nerve model. Saves
    the results to disk. This specific harmonic complex tone stimulus used in this simulation is from Guest and
    Oxenham (2021), although no acoustic noise is included in the stimulus. In this simulation, two separate 
    simulations are run --- one in which the even harmonics are simulated and another in which the odd harmonics
    are simulated. The results are then combined optimally across simulations.

    Args:
        model: model object from apcmodels.anf
        model_name (str): name of the model
        fs (int): sampling rate in Hz
    """
    # Define stimulus parameters
    F0s = 10**np.linspace(np.log10(280) - 0.2, np.log10(1400) + 0.1, 24)  # simulate 8th harmonic of F0s
    levels = [20, 40, 60]  # dB SPL, per component
    dur = 0.10  # seconds
    dur_ramp = 0.01  # seconds

    # Define model parameters
    cf_low = 5*F0s
    cf_high = 11*F0s
    n_cf = 40
    n_fiber_per_chan = round(((np.log10(11/5)/3)*18000)/n_cf)  # assume ~18k HSR fibers from 0.2 to 20 kHz

    # Encode parameters
    params = si.Parameters(dur=dur, dur_ramp=dur_ramp, fs=fs, n_cf=n_cf, delta_theta=[0.001], API=np.zeros(1),
                           n_fiber_per_chan=n_fiber_per_chan, model_name=model_name, comps=comps, fs_synapse=50e3)
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
    synth = ISOToneGuest2021_exp1a()
    stimuli = synth.synthesize_sequence(params)
    params.add_inputs(stimuli)

    # Construct simulation and run
    sim = model()
    results = sim.run(params, runfunc=dc.decode_ideal_observer(sim.simulate))

    # Compile results
    save_to_csv([res[0] for res in results], params, os.path.join('nofigure', 'parity', model_name + '_sim_dichotic_' + tag + '_AI.csv'), decoding_type='AI', model=model_name, roving_type='none', parity=tag)
    save_to_csv([res[1] for res in results], params, os.path.join('nofigure', 'parity', model_name + '_sim_dichotic_' + tag + '_RP.csv'), decoding_type='RP', model=model_name, roving_type='none', parity=tag)

# Loop through models and calculate FDLs for each model
simulate_f0dls_dichotic(anf.AuditoryNerveZilany2014, 'Zilany2014', 200e3, np.array([6, 7, 8, 9, 10]), 'both')
simulate_f0dls_dichotic(anf.AuditoryNerveZilany2014, 'Zilany2014', 200e3, np.array([6, 8, 10]), 'even')
simulate_f0dls_dichotic(anf.AuditoryNerveZilany2014, 'Zilany2014', 200e3, np.array([7, 9]), 'odd')


