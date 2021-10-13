"""
This script implements the simulations described in Figure 6 of Guest and Oxenham (2021).
"""
import apcmodels.simulation as si
import apcmodels.anf as anf
import apcmodels.decode as dc
import numpy as np
import os, sys
sys.path.append(os.getcwd())
from util.functions import ISOToneGuest2021, GEOMToneGuest2021, adjust_level
import matplotlib.pyplot as plt
plt.ion()


def get_ep_func(ratefunc):
    def inner(params):
        return [np.mean(ratefunc(x), axis=1) for x in params]
    return inner


def simulate(F0, model=anf.AuditoryNerveHeinz2001, model_name='Heinz2001', fs=500e3, delta=0.001, cf=2500, stim='iso', masker_interval=1):
   # Define stimulus parameters
    dur = 0.10  # seconds
    dur_ramp = 0.01  # seconds

    # Define model parameters
    n_fiber_per_chan = 1

    # Encode parameters
    params = si.Parameters(dur=dur, dur_ramp=dur_ramp, fs=fs, delta_theta=[delta], API=np.zeros(1),
                           n_fiber_per_chan=n_fiber_per_chan, model_name=model_name, F0=F0, level=30,
                           cf_low=cf, cf_high=cf, n_cf=1, F0_masker=F0*2**(masker_interval/12))

    # Adjust levels to be in dB re: threshold
    params.flatten()
    for ele in params:
        ele['nominal_level'] = ele['level']                                 # encode nominal level (dB re: threshold)
        ele['level'] = adjust_level(ele['F0']*8, ele['level'], model_name)  # encode actual level (dB SPL)

    # Encode increments
    params.increment({'F0': delta})  # increment F0

    # Synthesize stimuli
    if stim == 'iso':
        synth = ISOToneGuest2021()
    else:
        synth = GEOMToneGuest2021()
    stimuli = synth.synthesize_sequence(params)
    params.add_inputs(stimuli)

    # Construct simulation and run
    sim = model()
    return sim.run(params, runfunc=get_ep_func(sim.simulate), parallel=False)


# Show that GEOM vs ISO difference depends on delta size 
# Set parameters for figure
cf = 1975
resolution = 50

##### FIGURE 1
# Plot and sim
plt.figure()
deltas = 10**np.linspace(np.log10(1e-6), np.log10(5e1), resolution)
results = np.zeros_like(deltas)
for idx_delta, delta in enumerate(deltas):
    temp = simulate(280, delta=delta, cf=cf)
    results[idx_delta] = temp[0][1] - temp[0][0]

plt.plot(deltas, np.abs(results))
plt.xscale('log')
plt.yscale('log')

results = np.zeros_like(deltas)
for idx_delta, delta in enumerate(deltas):
    temp = simulate(280, delta=delta, stim='geom', cf=cf)
    results[idx_delta] = temp[0][1] - temp[0][0]

plt.plot(deltas, np.abs(results))
plt.xscale('log')
plt.yscale('log')
#plt.yscale('symlog')

##### FIGURE 2
# Plot and sim
plt.figure()
deltas = 10**np.linspace(np.log10(1e-6), np.log10(5e1), resolution)
results = np.zeros_like(deltas)
for idx_delta, delta in enumerate(deltas):
    temp = simulate(280, delta=delta, cf=cf)
    results[idx_delta] = temp[0][1] - temp[0][0]

plt.plot(deltas, results)
plt.xscale('log')
#plt.yscale('log')

results = np.zeros_like(deltas)
for idx_delta, delta in enumerate(deltas):
    temp = simulate(280, delta=delta, stim='geom', cf=cf)
    results[idx_delta] = temp[0][1] - temp[0][0]

plt.plot(deltas, results)
plt.xscale('log')
#plt.yscale('log')

##### FIGURE 3
# Plot and sim
plt.figure(figsize=(3,3))
deltas = 10**np.linspace(np.log10(1e-6), np.log10(5e1), resolution)
results = np.zeros_like(deltas)
for idx_delta, delta in enumerate(deltas):
    temp = simulate(280, delta=delta, cf=cf)
    results[idx_delta] = (temp[0][1] - temp[0][0])/delta

plt.plot(deltas, np.abs(results))
plt.xscale('log')
plt.yscale('log')

results = np.zeros_like(deltas)
for idx_delta, delta in enumerate(deltas):
    temp = simulate(280, delta=delta, stim='geom', cf=cf)
    results[idx_delta] = (temp[0][1] - temp[0][0])/delta

plt.plot(deltas, np.abs(results))
plt.xscale('log')
plt.yscale('log')
#plt.yscale('symlog')
plt.ylabel('Derivative estimate')
plt.xlabel('h')
plt.legend(['ISO', 'GEOM'])
plt.title('CF = ' + str(cf) + ' Hz')
plt.savefig(os.path.join('plots', 'supfig_maskers_derivative_estimate_figures.png'))