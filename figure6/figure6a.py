import apcmodels.simulation as si
import apcmodels.anf as anf
import numpy as np
from functions import ISOToneGuest2021, adjust_level
import matplotlib.pyplot as plt

# Configure figure
fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(10, 4))

# Run simulations
for idx_harmonic, harmonic in enumerate([8, 8.5]):
    for idx_F0, F0 in enumerate([280, 620, 1400]):
        for bool_noise, color, linewidth in zip([False, True], ['#fc8d62', 'gray'], [2, 1]):
            # Set up params
            params = si.Parameters(fs=int(200e3), fiber_type='hsr', n_cf=1, cf_low=F0*harmonic, cf_high=F0*harmonic,
                                   F0=F0, dur=(1/F0)*100, level=float(adjust_level(F0*8, 50, 'Zilany2014')),
                                   level_noise=float(adjust_level(F0*8, 40, 'Zilany2014')), ten=bool_noise)
            params.repeat(10)
            params.flatten()
            stim = ISOToneGuest2021()
            params.add_inputs(stim.synthesize_sequence(params))
            sim = anf.AuditoryNerveZilany2014()
            results = sim.run(params)
            means = np.squeeze(np.mean(results, axis=0))
            std = np.squeeze(np.std(results, axis=0))
            x_axis = np.linspace(0, means.shape[0]/200e3, means.shape[0])/(1/F0) - 35  # move onset to 35th cycle
            ax[idx_harmonic][idx_F0].plot(x_axis, means, color=color, linewidth=linewidth)
            ax[idx_harmonic][idx_F0].fill_between(x_axis, np.squeeze(means - std), np.squeeze(means + std),
                                   color=color, alpha=0.2, label='_nolegend_')
            # Handle axes
            if idx_harmonic == 0:
                ax[idx_harmonic][idx_F0].get_xaxis().set_visible(False)
            else:
                ax[idx_harmonic][idx_F0].set_xlabel('Time (periods of F0)')
            if idx_F0 > 0:
                ax[idx_harmonic][idx_F0].get_yaxis().set_visible(False)
            else:
                ax[idx_harmonic][idx_F0].set_ylabel('Firing rate (sp/s)')
            ax[idx_harmonic][idx_F0].set_xlim((0, 5))
            ax[idx_harmonic][idx_F0].set_ylim((0, 1000))
plt.tight_layout()
plt.savefig('plots/fig6.png')

res = []
for ii in range(100):
    x = ISOToneGuest2021().synthesize(F0=280, dur=1, dur_ramp=0.02, level=50, fs=int(48e3), ten=True, level_noise=50)
    X = 20*np.log10(np.abs(np.fft.fft(x)))
    res.append(X)
X = np.mean(np.array(res), axis=0)
plt.plot(X-np.max(X))
plt.xlim((0, 5000))
plt.ylim((-40, 5))