"""
The following functions are used throughout the scripts in this repo to synthesize the acoustic stimulus described in
Guest and Oxenham (2021) and to process some of the simulations.
"""
import apcmodels.synthesis as sy
import apcmodels.signal as sg
import numpy as np
from scipy.signal import sosfiltfilt, butter
from scipy.interpolate import interp1d


class ISOToneGuest2021(sy.Synthesizer):
    """
    Synthesizes the ISO stimulus in Guest and Oxenham (2021).
    """
    def __init__(self):
        super().__init__(stimulus_name='ISO Tone')

    def synthesize(self, F0, dur=0.350, dur_ramp=0.02, level=40, phase=0, fs=int(48e3), **kwargs):
        """
        Synthesizes a harmonic complex tone composed of all components of the F0 up to the 24000 Hz (the Nyquist
        frequency for the 48 kHz sampling rate used in the paper). Then, the tone is bandpass filtered from 5.5x to
        10.5x F0 using a zero-phase Butterworth filter.

        Arguments:
            F0 (float): F0 of the complex tone in Hz
            level (float, function): level per-component of the complex tone in dB SPL, can be either a float or
                a function that is evaluated for each component
            phase (float, function): phase offset applied to each component of the complex tone in degrees, can be
                either a float that is added to every component's phase or a function that is evaluated and whose
                output is then added to each component's phase
            dur (float): duration in seconds
            dur_ramp (float): duration of raised-cosine ramp in seconds
            fs (int): sampling rate in Hz

        Returns:
            output (array): complex tone stimulus
        """
        # Set up bandpass filter
        sos = butter(N=6, Wn=[F0 * 5.5 / (fs * 0.5), F0 * 10.5 / (fs * 0.5)], btype="band",
                     output="sos")
        # Create array of frequencies, levels, and phases
        freqs = np.arange(F0, 48000 / 2, F0)  # up to Nyquist for fs=48
        levels = np.zeros(len(freqs))
        for ii in range(len(levels)):
            if callable(level):
                levels[ii] = levels[ii] + level()
            else:
                levels[ii] = levels[ii] + level
        phases = np.zeros(len(freqs))
        for ii in range(len(phases)):
            if callable(phase):
                phases[ii] = phases[ii] + phase()
            else:
                phases[ii] = phases[ii] + phase
        # Synthesize, filter, and ramp complex tone signal
        signal = sg.complex_tone(freqs, levels, phases, dur, fs)
        signal = sosfiltfilt(sos, signal)
        signal = sg.cosine_ramp(signal, dur_ramp, fs)
        # Return
        return signal


def adjust_level(freq, level, model_name):
    """
    Accepts an input level in dB SPL and returns an adjusted level for the corresponding auditory nerve model.
    Adjustments are made based on estimates of absolute threshold for each nerve model made in
    nofigure/absolute_thresholds/.

    Parameters:
        freq (float): frequency at which the absolute threshold is estimated and used to adjust the level, in Hz
        level (float): input level in dB SPL
        model_name (str): either 'Heinz2001', 'Zilany2014', or 'Verhulst2018', indicates which model's absolute
            thresholds should be used in the adjustment
    """
    # Branch based on model
    if model_name == 'Heinz2001':
        cfs = np.array([200., 242.30553173, 293.55985352, 355.65588201,
                        430.88693801, 522.03144314, 632.45553203, 766.23736991,
                        928.31776672, 1124.68265038, 1362.58413812, 1650.80837054,
                        2000., 2423.05531726, 2935.59853524, 3556.55882008,
                        4308.86938006, 5220.31443137, 6324.55532034, 7662.37369911,
                        9283.17766723, 11246.82650381, 13625.84138116, 16508.08370536,
                        20000.])
        absolute_thresholds = np.array([11.75338957, 11.62867281, 11.68923837, 11.69510857, 11.71673396,
                                        11.78947423, 11.90428429, 12.0764798, 12.35406782, 12.76388996,
                                        13.21929098, 13.82237866, 14.72409746, 15.62127449, 16.58383025,
                                        17.33935012, 17.71101781, 17.87225817, 17.92201393, 17.9330261,
                                        17.93559906, 17.93581972, 17.93615349, 17.93697896, 17.93698821])
    else:
        return ValueError('Models other than Heinz (2001) are not yet implemented.')
    # Calculate correction and return
    adjustment = interp1d(np.log10(cfs), absolute_thresholds, kind='cubic')
    return level + adjustment(np.log10(freq))
