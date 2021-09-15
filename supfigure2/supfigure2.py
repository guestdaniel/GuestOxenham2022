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
from util.functions import adjust_level


class ISOToneGuest2021_exp1a_variable_harms(sy.Synthesizer):
    """ Synthesizes the ISO stimulus in Guest and Oxenham (2021), but with harmonic numbers that can be manipulated

    Simplified version of the stimulus in Guest and Oxenham (2021). This version is the version from Experiment 1a, and
    does not include acoustic masking noise.
    """
    def __init__(self):
        super().__init__(stimulus_name='ISO Tone')

    def synthesize(self, F0, dur=0.350, dur_ramp=0.02, level=None, phase=None, fs=int(48e3), h_low=6, **kwargs):
        """
        Synthesizes a harmonic complex tone composed of components 6-10 of the F0. This is the same stimulus as used
        in Experiment 1a.

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
        freqs = F0*(h_low+np.array([0, 1, 2, 3, 4]))
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


def simulate_supfigure2_f0dls_no_roving(model, model_name, fs):
    """
    Estimates F0 difference limens (F0DLs) using ideal observer analysis for a given auditory nerve model. Saves
    the results to disk. The simulations include phase randomization (specifically, each component from h_low-(h_low+4) has its
    starting phase randomized on each trial).

    Args:
        model: model object from apcmodels.anf
        model_name (str): name of the model
        fs (int): sampling rate in Hz
        n_rep (int): number of repetitions to perform for each simulation
    """
    # Define stimulus parameters
    F0s = 10**np.linspace(np.log10(280), np.log10(1400), 4)  # simulate 8th harmonic of F0s
    levels = [30]  # dB SPL
    dur = 0.10  # seconds
    dur_ramp = 0.01  # seconds
    h_lows = np.array([2, 6, 10])

    # Define model parameters
    cf_low = np.zeros((4, 3))
    cf_high = np.zeros((4, 3))
    for idx_F0, F0 in enumerate(F0s):
        for idx_h_low, h_low in enumerate(h_lows):
            cf_low[idx_F0, idx_h_low] = F0*(h_low-1)
            cf_high[idx_F0, idx_h_low] = F0*(h_low+5)
    n_cf = 40
    n_fiber_per_chan = round(((np.log10(11/5)/3)*18000)/n_cf)  # assume ~18k HSR fibers from 0.2 to 20 kHz

    # Construct delta_theta and API matrices/arrays
    delta_theta = [0.001]
    API = np.zeros((1, 1))

    # Encode parameters
    params = si.Parameters(dur=dur, dur_ramp=dur_ramp, fs=fs, n_cf=n_cf, delta_theta=delta_theta, API=API,
                           n_fiber_per_chan=n_fiber_per_chan, model_name=model_name)
    params.wiggle('F0', F0s)                                   # wiggle F0s
    params.wiggle('h_low', h_lows)
    params.stitch('cf_low', cf_low)                            # stitch cf_low (each F0 * h_low combo corresponds to a cf_low)
    params.stitch('cf_high', cf_high)                          # stitch cf_high (each F0 * h_low combo corresponds to a cf_high)
    params.wiggle('level', levels)                             # wiggle levels

    # Adjust levels to be in dB re: threshold
    params.flatten()
    for ele in params:
        ele['nominal_level'] = ele['level']                                 # encode nominal level (dB re: threshold)
        ele['level'] = adjust_level(ele['F0']*8, ele['level'], model_name)  # encode actual level (dB SPL)

    # Encode repeats and increments
    params.repeat(n_rep)
    params.increment({'F0': 0.001})                                 # increment F0

    # Synthesize stimuli
    synth = ISOToneGuest2021_exp1a_variable_harms()
    stimuli = synth.synthesize_sequence(params)
    params.add_inputs(stimuli)

    # Construct simulation and run
    sim = model()
    results = sim.run(params, runfunc=dc.decode_ideal_observer(sim.simulate))

    # Compile results
    save_to_csv([res[0] for res in results], params,
                'supfigure2/' + model_name + '_supfigure2_no_roving.csv', decoding_type='AI',
                model=model_name, roving_type='none')
    save_to_csv([res[1] for res in results], params,
                'supfigure2/' + model_name + '_supfigure2_no_roving_RP.csv', decoding_type='RP',
                model=model_name, roving_type='none')


def simulate_supfigure2_f0dls_phase_roving(model, model_name, fs, n_rep=10):
    """
    Estimates F0 difference limens (F0DLs) using ideal observer analysis for a given auditory nerve model. Saves
    the results to disk. The simulations include phase randomization (specifically, each component from h_low-(h_low+4) has its
    starting phase randomized on each trial).

    Args:
        model: model object from apcmodels.anf
        model_name (str): name of the model
        fs (int): sampling rate in Hz
        n_rep (int): number of repetitions to perform for each simulation
    """
    # Define stimulus parameters
    F0s = 10**np.linspace(np.log10(280), np.log10(1400), 4)  # simulate 8th harmonic of F0s
    levels = [30]  # dB SPL
    dur = 0.10  # seconds
    dur_ramp = 0.01  # seconds
    h_lows = np.array([2, 6, 10])

    # Define model parameters
    cf_low = np.zeros((4, 3))
    cf_high = np.zeros((4, 3))
    for idx_F0, F0 in enumerate(F0s):
        for idx_h_low, h_low in enumerate(h_lows):
            cf_low[idx_F0, idx_h_low] = F0*(h_low-1)
            cf_high[idx_F0, idx_h_low] = F0*(h_low+5)
    n_cf = 40
    n_fiber_per_chan = round(((np.log10(11/5)/3)*18000)/n_cf)  # assume ~18k HSR fibers from 0.2 to 20 kHz

    # Construct delta_theta and API matrices/arrays
    delta_theta = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001]
    API = np.zeros((6, 6))
    API[np.diag_indices(6)] = [0, 1/360**2, 1/360**2, 1/360**2, 1/360**2, 1/360**2]

    # Encode parameters
    params = si.Parameters(dur=dur, dur_ramp=dur_ramp, fs=fs, n_cf=n_cf, delta_theta=delta_theta, API=API,
                           n_fiber_per_chan=n_fiber_per_chan, model_name=model_name)
    params.wiggle('F0', F0s)                                   # wiggle F0s
    params.wiggle('h_low', h_lows)
    params.stitch('cf_low', cf_low)                            # stitch cf_low (each F0 * h_low combo corresponds to a cf_low)
    params.stitch('cf_high', cf_high)                          # stitch cf_high (each F0 * h_low combo corresponds to a cf_high)
    params.wiggle('level', levels)                             # wiggle levels

    # Adjust levels to be in dB re: threshold
    params.flatten()
    for ele in params:
        ele['nominal_level'] = ele['level']                                 # encode nominal level (dB re: threshold)
        ele['level'] = adjust_level(ele['F0']*8, ele['level'], model_name)  # encode actual level (dB SPL)

    # Encode repeats and increments
    params.repeat(n_rep)
    params.append('phase', lambda: np.random.uniform(0, 360, 5))  # append random phase
    params.increment({'F0': 0.001,                                  # increment F0
                      '(1)_phase': np.array([0.001, 0, 0, 0, 0]),   # increment phase of H6
                      '(2)_phase': np.array([0, 0.001, 0, 0, 0]),   # increment phase of H7
                      '(3)_phase': np.array([0, 0, 0.001, 0, 0]),   # increment phase of H8
                      '(4)_phase': np.array([0, 0, 0, 0.001, 0]),   # increment phase of H9
                      '(5)_phase': np.array([0, 0, 0, 0, 0.001])})  # increment phase of H10

    # Synthesize stimuli
    synth = ISOToneGuest2021_exp1a_variable_harms()
    stimuli = synth.synthesize_sequence(params)
    params.add_inputs(stimuli)

    # Construct simulation and run
    sim = model()
    results = sim.run(params, runfunc=dc.decode_ideal_observer(sim.simulate))

    # Compile results
    save_to_csv([res[0] for res in results], params,
                'supfigure2/' + model_name + '_supfigure2_phase_roving_AI.csv', decoding_type='AI',
                model=model_name, roving_type='phase')
    save_to_csv([res[1] for res in results], params,
                'supfigure2/' + model_name + '_supfigure2_phase_roving_RP.csv', decoding_type='RP',
                model=model_name, roving_type='phase')

# Loop through models and calculate FDLs for each model (for this figure, we only use Zilany et al. [2014])
for model, model_name, fs in zip([anf.AuditoryNerveZilany2014],
                                 ['Zilany2014'],
                                 [int(200e3)]):
    simulate_supfigure2_f0dls_no_roving(model, model_name, fs)
    simulate_supfigure2_f0dls_phase_roving(model, model_name, fs)
