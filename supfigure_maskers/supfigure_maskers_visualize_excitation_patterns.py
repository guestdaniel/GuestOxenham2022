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


def simulate_iso(F0, model, model_name, fs, delta=0.001, n_cf=40):
   # Define stimulus parameters
    dur = 0.10  # seconds
    dur_ramp = 0.01  # seconds

    # Define model parameters
    cf_low = 5*F0
    cf_high = 11*F0
    n_cf = n_cf
    n_fiber_per_chan = round(((np.log10(11/5)/3)*18000)/n_cf)  # assume ~18k HSR fibers from 0.2 to 20 kHz

    # Encode parameters
    params = si.Parameters(dur=dur, dur_ramp=dur_ramp, fs=fs, n_cf=n_cf, delta_theta=[delta], API=np.zeros(1),
                           n_fiber_per_chan=n_fiber_per_chan, model_name=model_name, F0=F0, level=30,
                           cf_low=cf_low, cf_high=cf_high)

    # Adjust levels to be in dB re: threshold
    params.flatten()
    for ele in params:
        ele['nominal_level'] = ele['level']                                 # encode nominal level (dB re: threshold)
        ele['level'] = adjust_level(ele['F0']*8, ele['level'], model_name)  # encode actual level (dB SPL)

    # Encode increments
    params.increment({'F0': delta})  # increment F0

    # Synthesize stimuli
    synth = ISOToneGuest2021()
    stimuli = synth.synthesize_sequence(params)
    params.add_inputs(stimuli)

    # Construct simulation and run
    sim = model()
    return sim.run(params, runfunc=dc.decode_ideal_observer(sim.simulate)), sim.run(params, runfunc=get_ep_func(sim.simulate)), params


def simulate_geom(F0, model, model_name, fs, n_rep=10, masker_interval=1, delta=0.001, n_cf=40):
   # Define stimulus parameters
    dur = 0.10  # seconds
    dur_ramp = 0.01  # seconds

    # Define model parameters
    cf_low = 5*F0
    cf_high = 11*F0
    n_cf = n_cf
    n_fiber_per_chan = round(((np.log10(11/5)/3)*18000)/n_cf)  # assume ~18k HSR fibers from 0.2 to 20 kHz

    # Encode parameters
    params = si.Parameters(dur=dur, dur_ramp=dur_ramp, fs=fs, n_cf=n_cf, delta_theta=[delta], API=np.zeros(1),
                           n_fiber_per_chan=n_fiber_per_chan, model_name=model_name, F0=F0, F0_masker=F0*2**(masker_interval/12), 
                           level=30, cf_low=cf_low, cf_high=cf_high)

    # Adjust levels to be in dB re: threshold
    params.flatten()
    for ele in params:
        ele['nominal_level'] = ele['level']                                 # encode nominal level (dB re: threshold)
        ele['level'] = adjust_level(ele['F0']*8, ele['level'], model_name)  # encode actual level (dB SPL)

    # Encode increments
    params.increment({'F0': delta})  # increment F0

    # Synthesize stimuli
    synth = GEOMToneGuest2021()
    stimuli = synth.synthesize_sequence(params)
    params.add_inputs(stimuli)

    # Construct simulation and run
    sim = model()
    return sim.run(params, runfunc=dc.decode_ideal_observer(sim.simulate)), sim.run(params, runfunc=get_ep_func(sim.simulate)), params

# Set parameters
n_cf = 100
cfs = 10**np.linspace(np.log10(280*5), np.log10(280*11), n_cf)

# Show that GEOM vs ISO difference depends on delta size 
deltas = [1e-2, 1e-1, 1e0, 1e1]
fig, axs = plt.subplots(1, 4, figsize=(8, 3), sharey=True)
for idx_delta, delta in enumerate(deltas):
    # Draw horizontal line at 0
    axs[idx_delta].plot([1400, 3100], [0, 0], color='black')
    # Draw ISO excitation pattern
    thresholds_iso, eps_iso, params_iso = simulate_iso(280, anf.AuditoryNerveHeinz2001, 'Heinz2001', 500e3, delta=delta, n_cf=n_cf)
    ep_diff = (eps_iso[0][1]-eps_iso[0][0])/delta
    axs[idx_delta].plot(cfs, ep_diff)

    # Draw GEOM excitation pattern
    thresholds_geom, eps_geom, params_geom = simulate_geom(280, anf.AuditoryNerveHeinz2001, 'Heinz2001', 500e3, delta=delta, n_cf=n_cf)
    ep_diff = (eps_geom[0][1]-eps_geom[0][0])/delta
    axs[idx_delta].plot(cfs, ep_diff)

    # Draw harmonic indicators
    for harm in [6, 7, 8, 9, 10]:
        axs[idx_delta].plot([280*harm, 280*harm], [-10, 10], color='gray', linestyle='dashed')
        axs[idx_delta].arrow(x=harm*280, y=5, dx=30, dy=0)

    # Set ylims
    axs[idx_delta].set_ylim((-10, 10))

    # Set titles
    axs[idx_delta].set_title('h = ' + str(delta) + ' Hz')
    axs[idx_delta].set_xlabel('CF (Hz)')

axs[0].set_ylabel('Derivative estimate ([sp/s]/Hz)')
plt.savefig(os.path.join('plots', 'excitation_pattern_differences.png'), dpi=300)

# Show instability in GEOM finite difference estimates
#deltas = [1e-2, 1e-1, 1e0, 1e1]
#fig, axs = plt.subplots(1, 2)
#for idx_delta, delta in enumerate(deltas):
#    thresholds_iso, eps_iso, params_iso = simulate_iso(280, anf.AuditoryNerveHeinz2001, 'Heinz2001', 500e3, delta=delta, n_cf=n_cf)
#    ep_diff = eps_iso[0][1]-eps_iso[0][0]
#    axs[0].plot(cfs, ep_diff/np.max(np.abs(ep_diff)))

#    thresholds_geom, eps_geom, params_geom = simulate_geom(280, anf.AuditoryNerveHeinz2001, 'Heinz2001', 500e3, delta=delta, n_cf=n_cf)
#    ep_diff = eps_geom[0][1]-eps_geom[0][0]
#    axs[1].plot(cfs, ep_diff/np.max(np.abs(ep_diff)))
#axs[1].legend(deltas)