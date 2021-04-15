"""
This script implements Figure 0 of Guest and Oxenham (2021).
"""
import apcmodels.synthesis as sy
import apcmodels.simulation as si
import apcmodels.anf as anf
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os, sys
sys.path.append(os.getcwd())


def sim_and_plot_neurogram(f0, xlow=20, xhigh=25):
    """ Plots a simulated auditory nerve response for a harmonic complex tone.

    The complex tone is synthesized in sine phase at 50 dB SPL per-component with components 4-13.

    Args:
        f0 (float): Fundamental frequency of the HCT, in Hz.
        xlow (int): the lower limit of the x-axis, in periods of the tone
        xhigh (int); the upper limit of the x-axis, in periods of the tone
    """
    # Simulate neurogram
    params = si.Parameters(f0=f0, level=50, fs=int(200e3), cf_low=f0*0.5, cf_high=f0*15, h_low=4, h_high=13, n_cf=200,
                           dur=0.10, dur_ramp=0.01, fiber_type='msr')
    stimulus = sy.ComplexTone().synthesize_sequence(params)
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
    fig, axs = plt.subplots(3, 3, gridspec_kw={'width_ratios': [1, 4, 1], 'height_ratios': [1, 4, 1]})
    # Null out corners
    axs[0, 0].axis('off')
    axs[2, 2].axis('off')
    axs[0, 2].axis('off')
    axs[2, 0].axis('off')

    # Plot surface (middle)
    axs[1, 1].pcolormesh(t, f, resp[0], shading='auto')
    axs[1, 1].get_xaxis().set_visible(False)
    axs[1, 1].get_yaxis().set_visible(False)
    axs[1, 1].set_xlim((xlow, xhigh))

    # Plot acoustic stimulus (top)
    axs[0, 1].plot(t, stim, color='black')
    axs[0, 1].get_xaxis().set_visible(False)
    axs[0, 1].set_xlim((xlow, xhigh))
    axs[0, 1].set_ylabel('Amplitude (Pa)', rotation='horizontal', ha='right')

    # Plot average neural response (bottom)
    idx = np.argmin(np.abs(f - f0*4))
    h4 = axs[2, 1].plot(t, resp[0][idx, :], color='slateblue')
    idx = np.argmin(np.abs(f - f0*12))
    h12 = axs[2, 1].plot(t, resp[0][idx, :], color='cornflowerblue')
    axs[2, 1].set_xlim((xlow, xhigh))
    axs[2, 1].set_xlabel('Normalized time (periods)')
    axs[2, 1].legend(['H4', 'H12'], framealpha=1)
    axs[2, 1].set_ylabel('Firing rate', rotation='horizontal')
    axs[2, 1].yaxis.tick_right()
    axs[2, 1].yaxis.set_label_position('right')
    axs[2, 1].yaxis.set_label_coords(1.20, 0.40)

    # Plot spectrum (left)
    X = 20*np.log10(np.abs(np.fft.fft(stim, n=stim.shape[0]*10)))
    X = X - np.max(X)
    faxis = np.linspace(0, params[0]['fs'], X.shape[0])
    rect = patches.Rectangle((-10, f0*5.5), 12, f0*5, edgecolor='none', facecolor='gold')
    axs[1, 0].add_patch(rect)
    axs[1, 0].plot(X, faxis, color='black')
    axs[1, 0].plot([-10, 0], [f0, f0], 'r--')
    axs[1, 0].set_ylim(params[0]['cf_low'], params[0]['cf_high'])
    axs[1, 0].set_xlim((-10, 5))
    axs[1, 0].set_ylabel('Frequency or CF (Hz)')
    axs[1, 0].set_xlabel('Level (dB)')

    # Plot excitation pattern (right)
    axs[1, 2].plot(np.mean(resp[0], axis=1), f)
    idx = np.argmin(np.abs(f - f0*4))
    axs[1, 2].plot([0, np.mean(resp[0], axis=1)[idx]], [f0*4, f0*4], color='slateblue', linestyle='dashed')
    idx = np.argmin(np.abs(f - f0*12))
    axs[1, 2].plot([0, np.mean(resp[0], axis=1)[idx]], [f0*12, f0*12], color='cornflowerblue', linestyle='dashed')
    axs[1, 2].set_ylim(params[0]['cf_low'], params[0]['cf_high'])
    axs[1, 2].get_yaxis().set_visible(False)
    axs[1, 2].set_xlabel('Firing rate')
    axs[1, 2].xaxis.tick_top()
    axs[1, 2].xaxis.set_label_position('top')

    # Tighten layout
    plt.tight_layout(pad=0.1, h_pad=-1, w_pad=-10)

# Plot F0 = 280 Hz
sim_and_plot_neurogram(280, 20, 25)
plt.savefig('plots/fig0a.png')

# Plot F0 = 1400 Hz
sim_and_plot_neurogram(1400, 20*4, 20*4+5)
plt.savefig('plots/fig0b.png')
