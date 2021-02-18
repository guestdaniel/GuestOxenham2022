"""
This script estimates Q10 based on clicks
"""
import apcmodels.synthesis as sy
import apcmodels.simulation as si
import apcmodels.anf as anf
import apcmodels.signal as sg
import numpy as np
import config as cfg
import os
import sys
sys.path.append('/home/daniel/apc_code/scripts/Verhulstetal2018Model')
from run_model2018 import Verhulst2018Cochlea


def calc_erb(freq):
    """
    Calculates the ERB at a given frequency

    Arguments:
        freq (float): frequency in Hz

    Returns:
        ERB (float): ERB in Hz
    """
    freq = freq/1000  # transform from Hz to kHz
    return 24.7 * (4.37 * freq + 1)  # equation from Moore and Glasberg (1983)


def synthesize_click(dur_click=0.08 / 1000, dur_pre_click=0.01, dur_post_click=0.06, fs=int(200e3)):
    """
    Synthesizes a single click embedded in silence

    Arguments:
        dur_click (float): duration of click in seconds
        dur_stim (float): duration of stimulus in seconds. The click is embedded in the middle of the stimulus
        fs (int): sampling rate in Hz

    Returns:
        output (array): pure tone
    """
    dur_stim = dur_pre_click + dur_click + dur_post_click
    # Create time axis
    t = np.linspace(0, dur_stim, int(dur_stim * fs))
    # Create click
    click = np.zeros(shape=t.shape)
    click[np.logical_and(t > dur_pre_click, t < (dur_pre_click+dur_click))] = 1
    return click


def calculate_verhulst2018_bm_response(_input, fs, **kwargs):
    """
    Implements Verhulst, Altoe, and Vasilikov (2018) peripheral model and returns basilar membrane response

    Arguments:
        _input (ndarray): 1-dimensional ndarray containing an acoustic stimulus in pascals
        fs (int): sampling rate in Hz

    Returns:
        output (ndarray): output array of basilar membrane response, of shape (n_cf, n_sample)

    Warnings:
        - Note that arguments passed to **kwargs are discarded silently
    """
    # Run firing rate simulations
    vm, fs_res, cfs_model = Verhulst2018Cochlea(_input, fs)
    vm = np.flip(vm, axis=1)  # flip tonotopic axis left-right
    vm = vm.T  # transpose to (n_cf, n_sample)

    return vm, fs_res, np.flip(cfs_model)


def QERB_calculation(bmm, cf, fs):
    samples = bmm.shape[0]
    half = samples / 2
    ener = 0
    F = (2 * abs(np.fft.fft(bmm)) / samples) ** 2
    max_val = F.max(0)
    for j in range(int(half) + 1):
        ener = ener + F[j]
    BW = (ener / max_val) * fs / samples
    Q = cf / BW
    return Q


def QERB_calculation(bmm, cfs, fs):
    central = cfs.shape[0]
    samples = bmm.shape[1]
    half = samples / 2
    F = np.zeros((samples, central))
    G = np.zeros((samples, central))
    max_val = np.zeros(central)
    ener = np.zeros(central)
    BW = np.zeros(central)
    QdB = np.zeros(central)

    for i in range(int(central)):
        F[:, i] = (2 * abs(np.fft.fft(bmm[i, :])) / samples) ** 2
        max_val[i] = F.max(0)[i]
        for j in range(int(half) + 1):
            ener[i] = ener[i] + F[j, i]
        # ener[i] = (F.sum(0)[i])/2
        BW[i] = (ener[i] / max_val[i]) * fs / samples
        QdB[i] = cfs[i] / BW[i]
    return QdB


# Simulate Verhulst basilar membrane response to 30 dB SPL click
fs = int(200e3)
click = synthesize_click(fs=fs)
click = sg.scale_dbspl(click, 50)
vm, _, model_cfs = calculate_verhulst2018_bm_response(click, fs=int(200e3))

# Calculate energy in impulse response at each CF we care about
cfs = 10**np.linspace(np.log10(200), np.log10(20000), 25)  # CFs for which we will measure q10
q10s = []
#for cf in cfs:
    # Extract the impulse response at this CF and calculate its power
    #bm_impulse_response = vm[np.argmin(np.abs(model_cfs-cf)), :]
    #power_in_impulse_response = np.sum(bm_impulse_response**2)  # calculate power as int x^2(t) dt, Parseval's theorem
    #q10s.append(QERB_calculation(bm_impulse_response, cf, fs))
    # Calculate the power that would enter into an ERB at this CF
    #X = np.fft.fft(click)  # calculate the spectrum of the click
    #f = np.linspace(0, fs, len(click))  # calculate the frequency axis for X
    #erb = calc_erb(cf)  # calculate ERB at this cf
    #X[np.logical_or(np.logical_not(np.logical_and(f >= cf-erb/2, f <= cf+erb/2)),  # zero out spectrum not in passband
    #                np.logical_not(np.logical_and(f >= (fs-cf)-erb/2, f <= (fs-cf)+erb/2)))] = 0
    #power_in_erb = np.sum(np.abs(X)**2)

#cf_test = np.array([1000, 2000, 4000, 8000])
plt.plot(model_cfs, QERB_calculation(vm, model_cfs, fs)/2)