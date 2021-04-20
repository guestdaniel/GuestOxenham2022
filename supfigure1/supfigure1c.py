import apcmodels.simulation as si
import apcmodels.anf as anf
import numpy as np
from util.functions import ComplexToneCedolin2005
import matplotlib.pyplot as plt
import matplotlib
from scipy.interpolate import make_interp_spline
import os, sys
sys.path.append(os.getcwd())


def plot_excitation_pattern_Cedolin_2005(F0s, cf, ax, first):
    """
    Simulates mean responses for a single auditory nerve fiber tuned to a particular cf as a function of the F0 of
    the stimulus. Simulations are performed at multiple levels (20, 30, and 40) dB SPL.

    Parameters:
        F0s (np.ndarray): array of F0 values to test, of sample (n_F0, )
        cf (float): characteristic frequency value to test
        ax (axis): axis object to display on
        first (bool): if first, we don't label x-axis
    """
    # Choose parameters to run
    levels = [20, 30, 40]  # per-component level in dB SPL

    # Setup params
    params = si.Parameters(fs=int(200e3), fiber_type='hsr', n_cf=1, cf_low=cf, cf_high=cf, species='cat')
    params.wiggle('level', levels)
    params.wiggle('F0', F0s)  # params is now shape of (3, 50)

    # Synthesize stimuli and add to params
    stim = ComplexToneCedolin2005()
    params.add_inputs(stim.synthesize_sequence(params))

    # Select model and run
    sim = anf.AuditoryNerveZilany2014()
    results = sim.run(params, runfunc=lambda x: np.mean(sim.simulate(x)))  # wrap default rate simulation in np.mean

    # Plot
    for ii in range(3):
        # Plot underlying data with scatterplot
        ax.scatter(cf / F0s, results[ii, :], color=matplotlib.cm.get_cmap('Blues')((ii + 1) * 0.25))
        # Interpolate data at higher resolution
        spl = make_interp_spline(F0s, results[ii, :], k=5)
        F0s_new = 10 ** np.linspace(np.log10(np.min(F0s)), np.log10(np.max(F0s)), num=500)
        results_smooth = spl(F0s_new)
        ax.plot(cf / F0s_new, results_smooth, color=matplotlib.cm.get_cmap('Blues')((ii + 1) * 0.25))

    # Add labels
    ax.set_title(str(cf) + ' Hz')
    if first:
        ax.get_xaxis().set_visible(False)
        ax.legend(['20', '30', '40'], title='Level (dB SPL)')
    else:
        ax.set_xlabel('Neural harmonic number (CF/F0)')
    ax.set_ylabel('Firing rate (sp/s)')

# Plot supfigure1c for both cf=952 and cf=4026
fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(4, 4))
plot_excitation_pattern_Cedolin_2005(10**np.linspace(np.log10(106), np.log10(474), 50), 952, ax[0], True)
plot_excitation_pattern_Cedolin_2005(10**np.linspace(np.log10(447), np.log10(2012), 50), 4026, ax[1], False)
plt.tight_layout()
# Save plot to disk
plt.savefig('plots/supfig1c.png')