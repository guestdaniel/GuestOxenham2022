import apcmodels.synthesis as sy
import apcmodels.simulation as si
import apcmodels.anf as anf
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os, sys
sys.path.append(os.getcwd())
import matplotlib
matplotlib.use('Agg')
from util.functions import ISOToneGuest2021, GEOMToneGuest2021, DBLToneGuest2021


def sim_and_plot_neurogram(f0, stimulus, xlow=20, xhigh=25):
    """ Plots a simulated auditory nerve response for a harmonic complex tone.

    The complex tone is synthesized in sine phase at 50 dB SPL per-component with components 4-13.

    Args:
        f0 (float): Fundamental frequency of the HCT, in Hz.
        xlow (int): the lower limit of the x-axis, in periods of the tone
        xhigh (int); the upper limit of the x-axis, in periods of the tone
    """
    # Simulate neurogram
    if stimulus == ISOToneGuest2021:
        params = si.Parameters(F0=f0, level=50, fs=int(200e3), cf_low=f0*0.5, cf_high=f0*15, n_cf=200,
                            dur=0.10, dur_ramp=0.01, fiber_type='msr')
    elif stimulus == GEOMToneGuest2021:
        params = si.Parameters(F0=f0, F0_masker=f0*2**(0.5/12), level=50, level_masker=50, fs=int(200e3), cf_low=f0*0.5, cf_high=f0*15, n_cf=200,
                            dur=0.10, dur_ramp=0.01, fiber_type='msr')
    else:
        params = si.Parameters(F0=f0, F0_masker_1=f0*2**(-6/12), F0_masker_2=f0*2**(8/12), level=50, level_masker_1=50, level_masker_2=50, fs=int(200e3), cf_low=f0*0.5, cf_high=f0*15, n_cf=200,
                            dur=0.10, dur_ramp=0.01, fiber_type='msr')
    stimulus = stimulus().synthesize_sequence(params)
    params.add_inputs(stimulus)
    sim = anf.AuditoryNerveZilany2014()
    resp = sim.run(params)

    # Calculate various axes
    t = np.linspace(0, resp[0].shape[1]/params[0]['fs'], resp[0].shape[1])
    t = t / (1/f0) - xlow  # adjust by requested xlow so that the x-axis will start at 0
    xhigh = xhigh - xlow
    xlow = xlow - xlow
    f = 10**np.linspace(np.log10(params[0]['cf_low']), np.log10(params[0]['cf_high']), params[0]['n_cf'])

    # Get stimulus and add zero padding to it
    stim = np.concatenate([np.zeros(int(params[0]['fs']*0.005)),
                                        params[0]['_input'],
                                        np.zeros(int(params[0]['fs']*0.040))])

    # Generate plot
    fig, axs = plt.subplots(2, 3, gridspec_kw={'width_ratios': [1, 4, 1], 'height_ratios': [1, 4]}, figsize=(7, 5))
    # Null out corners
    axs[0, 0].axis('off')
    axs[0, 2].axis('off')

    # Plot surface (middle)
    axs[1, 1].pcolormesh(t, f, resp[0], shading='auto')
    #axs[1, 1].get_xaxis().set_visible(False)
    axs[1, 1].set_xlabel("Normalized time (periods)")
    axs[1, 1].get_yaxis().set_visible(False)
    axs[1, 1].set_xlim((xlow, xhigh))

    # Plot acoustic stimulus (top)
    axs[0, 1].plot(t, stim, color='black')
    axs[0, 1].get_xaxis().set_visible(False)
    axs[0, 1].set_xlim((xlow, xhigh))
    axs[0, 1].set_ylabel('Amplitude (Pa)', rotation='horizontal', ha='right')

    # Plot spectrum (left)
    X = 20*np.log10(np.abs(np.fft.fft(stim, n=stim.shape[0]*10)))
    X = X - np.max(X)
    faxis = np.linspace(0, params[0]['fs'], X.shape[0])
    axs[1, 0].plot(X, faxis, color='black')
    axs[1, 0].plot([-10, 0], [f0, f0], 'r--')
    axs[1, 0].set_ylim(params[0]['cf_low'], params[0]['cf_high'])
    axs[1, 0].set_xlim((-10, 5))
    axs[1, 0].set_ylabel('Frequency or CF (Hz)')
    axs[1, 0].set_xlabel('Level (dB)')

    # Plot excitation pattern (right)
    axs[1, 2].plot(np.mean(resp[0], axis=1), f)
    idx = np.argmin(np.abs(f - f0*4))
    idx = np.argmin(np.abs(f - f0*12))
    idx = np.argmin(np.abs(f - f0*5.5))
    idx = np.argmin(np.abs(f - f0*9.5))
    axs[1, 2].set_ylim(params[0]['cf_low'], params[0]['cf_high'])
    axs[1, 2].get_yaxis().set_visible(False)
    axs[1, 2].set_xlabel('Firing rate')
    axs[1, 2].xaxis.tick_top()
    axs[1, 2].xaxis.set_label_position('top')

    # Tighten layout
    plt.tight_layout(pad=0.1, h_pad=-1, w_pad=-10)

# Plot F0 = 280 Hz
sim_and_plot_neurogram(280, ISOToneGuest2021, 20, 25)
plt.savefig('plots/s5_text_figa1.png')
sim_and_plot_neurogram(1400, ISOToneGuest2021, 20*4, 20*4+5)
plt.savefig('plots/s5_text_figa2.png')
sim_and_plot_neurogram(280, GEOMToneGuest2021, 20, 25)
plt.savefig('plots/s5_text_figa3.png')
sim_and_plot_neurogram(1400, GEOMToneGuest2021, 20*4, 20*4+5)
plt.savefig('plots/s5_text_figa4.png')
sim_and_plot_neurogram(280, DBLToneGuest2021, 20, 25)
plt.savefig('plots/s5_text_figa5.png')
sim_and_plot_neurogram(1400, DBLToneGuest2021, 20*4, 20*4+5)
plt.savefig('plots/s5_text_figa6.png')


