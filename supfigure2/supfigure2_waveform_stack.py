"""
This script implements the simulations described in Figure 6b of Guest and Oxenham (2021).
"""
import apcmodels.simulation as si
import apcmodels.anf as anf
import apcmodels.synthesis as sy
import apcmodels.decode as dc
import apcmodels.signal as sg
from apcmodels.util import save_to_csv
import numpy as np
import os, sys
sys.path.append(os.getcwd())
from util.functions import adjust_level
import matplotlib.pyplot as plt


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


# Now, look at responses off-CF as function of relative phase
def sim_one_response2(F0, h_low, n_harm, cf=None, n_samp=50):
    params = si.Parameters(dur=0.10, dur_ramp=0.01, level=30, fs=200e3, n_cf=1, F0=F0,
                           h_low=h_low, n_harm=n_harm,
                           n_fiber_per_chan=50, cf_low=cf, cf_high=cf)

    for ele in params:
        # encode nominal level (dB re: threshold)
        ele['nominal_level'] = ele['level']
        # encode actual level (dB SPL)
        ele['level'] = adjust_level(ele['F0']*8, ele['level'], 'Zilany2014')

    phase = [np.random.uniform(0, 360, n_harm) for samp in range(n_samp)]
    params.wiggle('phase', phase)

    # Synthesize stimuli
    synth = ISOToneGuest2021_exp1a_variable_harms()
    stimuli = synth.synthesize_sequence(params)
    params.add_inputs(stimuli)

    # Construct simulation and run
    sim = anf.AuditoryNerveZilany2014()
    results = sim.run(params)

    # Return
    return(results)


# Write fancier function to plot good version of these plots
def plot_fancy_stack(F0, n_samp=50):
    h_lows = [2, 6, 10]
    fig, axs = plt.subplots(3, 1, sharex=True, sharey=True)
    for ax, h_low in zip(axs, h_lows):
        results = sim_one_response2(F0, h_low, 5, F0*(h_low+0.5), n_samp)
        mean_response = np.squeeze(np.mean(results).T)
        sd_response = np.squeeze(np.std(results).T)
        t = np.linspace(0, len(mean_response)/200e3, len(mean_response))
        ax.plot(t, mean_response)
        ax.fill_between(t, mean_response - sd_response, mean_response + sd_response,
                        alpha=0.2, label='_nolegend_')
        ax.set_xlim((0.03, 0.04))
        # Clean up axis
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
    plt.ylabel('Firing rate (sp/s)')
    plt.xlabel('Time (s)')
    plt.show()

plot_fancy_stack(1400, n_samp=1000)