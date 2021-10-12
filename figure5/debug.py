"""
This script implements the simulations described in Figure 5a of Guest and Oxenham (2021).
"""
import apcmodels.simulation as si
import apcmodels.anf as anf
import numpy as np
import itertools
import os, sys
sys.path.append(os.getcwd())
from util.functions import DBLToneGuest2021
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import util as cfg
from util.functions import adjust_level


def resp(f0, tmr, level_noise, cf):
    """ Synthesizes complex tone mixture and simulates interspike interval histogram response.

    Args:
        f0 (float): cf of the complex tone
        neural_harm_nums (np.ndarray): array of neural harmonic numbers (CF/F0) to test. This is performed by keeping
            the F0 of the tone fixed while manipulating the range of CFs tested.
        tmr (float): Target-to-masker ratio of the target complex tone vs the masker complex tone
        n_repeat (int): number of repeats to run for each CF/F0 to test

    Returns:
        histograms (list): list of histograms counts for ISI histogram for each neural harmonic number
    """
    # Setup params [simulation]
    params = si.Parameters(fs=int(400e3), fs_synapse=200e3, n_cf=1, anf_num=(1, 0, 0),                        # model
                           F0=f0, F0_masker_1=f0*2**(-5.5/12), F0_masker_2=f0*2**(6/12),                      # F0s
                           level=50, level_masker_1=50-tmr, level_masker_2=50-tmr,                            # levels
                           phase=0, phase_masker_1=0, phase_masker_2=0,                                       # phases
                           ten=True, level_noise=level_noise, cf_low=cf, cf_high=cf)                                                          # noise

    # Synthesize stimuli and add to params
    stim = DBLToneGuest2021()
    params.add_inputs(stim.synthesize_sequence(params))

    # Select model and run
    sim = anf.AuditoryNerveZilany2014()
    results = sim.run(params)

    # Return
    return results


plt.figure()
for level_noise in [-100, 37]:
    x = resp(1400, 50, level_noise, 1400*7.5)
    plt.subplot(3, 1, 1)
    plt.plot(np.squeeze(x[0]))
    plt.subplot(3, 1, 2)
    f = np.linspace(0, 400e3, len(np.squeeze(x[0])))
    X = 20*np.log10(np.abs(np.fft.fft(np.squeeze(x[0]))))
    plt.plot(f, X-np.max(X))
    plt.xlim((0, 20e3))
    for harm in np.arange(1, 21):
        plt.arrow(harm*1400, 10, 0, -5, color='red')
    plt.ylim((-50, 15))
    plt.subplot(3, 1, 3)
    y = np.squeeze(x[0])
    lag = np.linspace(0, len(y)/400e3, len(y))
    plt.plot(lag*1400, np.fft.ifft(np.abs(np.fft.fft(y))**2))
    plt.xlim((0, 10))
