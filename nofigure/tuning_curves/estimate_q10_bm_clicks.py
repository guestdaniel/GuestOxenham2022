"""
This script estimates peripheral tuning for the Verhulst et al. (2018) peripheral model. This is done by presenting
low-level clicks to the basilar membrane stage of the model and then computing QERB from the energy under the power
spectrum of the impulse response at each CF. These values are then converted to Q10 using the conversion formula
suggested Verschooten et al. that Q10 = QERB/1.83.
"""

import os, sys
sys.path.append(os.getcwd())
import util as cfg
import numpy as np
from apcmodels.external.verhulst2018.run_verhulst import run_verhulst2018_cochlea


def synthesize_click(dur_click=0.08 / 1000, dur_pre_click=1/1000, dur_post_click=30/1000, fs=int(200e3)):
    """
    Synthesizes a single click embedded in silence

    Arguments:
        dur_click (float): duration of click in seconds
        dur_pre_click (float): duration of the time before the click in seconds
        dur_post_click (float): duration of the time after the click in seconds
        fs (int): sampling rate in Hz

    Returns:
        output (array): pure tone
    """
    # Calculate the duration of the stimulus
    dur_stim = dur_pre_click + dur_click + dur_post_click
    # Create time axis
    t = np.linspace(0, dur_stim, int(dur_stim * fs))
    # Create click
    click = np.zeros(shape=t.shape)
    click[np.logical_and(t > dur_pre_click, t < (dur_pre_click+dur_click))] = 1
    return click


def calculate_verhulst2018_bm_response(_input, fs, **kwargs):
    """
    Returns basilar membrane motion from Verhulst, Altoe, and Vasilikov (2018) peripheral model

    Arguments:
        _input (ndarray): 1-dimensional ndarray containing an acoustic stimulus in pascals
        fs (int): sampling rate in Hz

    Returns:
        output (ndarray): output array of basilar membrane response, of shape (n_cf, n_sample)

    Warnings:
        - Note that arguments passed to **kwargs are discarded silently
        - fs is not currently used
    """
    # Run firing rate simulations
    vm, cfs = run_verhulst2018_cochlea(_input, fs)
    vm = np.flip(vm, axis=1)  # flip tonotopic axis left-right
    vm = vm.T  # transpose to (n_cf, n_sample)

    return vm, np.flip(cfs)


def calc_Q(vm, cfs, fs):
    """
    Calculates QERB at each channel in vm under the assumption that vm reflects basilar membrane impulse responses to a
    short click stimulus.

    Arguments:
        vm (ndarray): array of basilar membrane vibration, of shape (n_channel, n_sample)
        cfs (ndarray): array of characteristic frequencies for each channel
        fs (int): sampling rate in Hz

    Returns:
        Q (list): QERB at each CF
    """
    Q = list()  # create empty list to store Q
    # Loop through channels
    for ii in range(vm.shape[0]):
        X = (np.abs(np.fft.fft(vm[ii, :])))**2  # calculate power spectrum at CF
        X = X[0:int(len(X)/2)]  # take LHS only of spectrum
        ERB = np.sum(X)/np.max(X) * fs / vm.shape[1]  # Calculate ERB
        # Note: note that sum of energy divided by max energy is ERB
        Q.append(cfs[ii]/ERB)
    return Q


# Simulate Verhulst basilar membrane response to 30 dB SPL click
fs = int(500e3)
click = synthesize_click(fs=fs)
p0 = 2e-5
L = 40
click = click*2*np.sqrt(2)*p0*10**(L/20)
vm, model_cfs = calculate_verhulst2018_bm_response(click, fs)

# Calculate QERB values from the basilar membrane responses
q_erb = calc_Q(vm, model_cfs, fs)

# Convert the QERB values to Q10 values
q_10 = np.array(q_erb)/1.83

# Subset the q_10 values to only include those that were tested in the other models
cfs = np.load('nofigure/tuning_curves/cfs.npy')  # load the tested CF values
idxs = [np.argmin(np.abs(model_cfs-cf)) for cf in cfs]  # locate the indices for Verhulst CFs closest to tested CFs
q_10 = q_10[idxs]  # subset the q_10 values

# Save these values out to disk
np.save('nofigure/tuning_curves/Verhulst2018_q10s' + '.npy', q_10)