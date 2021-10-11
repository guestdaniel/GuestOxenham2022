import apcmodels.synthesis as sy
import apcmodels.simulation as si
import apcmodels.anf as anf
import apcmodels.decode as dc
from apcmodels.util import save_to_csv
import numpy as np
import os, sys
sys.path.append(os.getcwd())
from util.functions import ISOToneGuest2021_exp1a, adjust_level
import util as cfg
import matplotlib.pyplot as plt


def get_information_curves(ratefunc):

    def inner(params):
        # Pull parameters from encoded list/dict of parameters
        fs = dc.find_parameter(params, 'fs')
        delta_theta = dc.find_parameter(params, 'delta_theta')
        n_fiber_per_chan = dc.find_parameter(params, 'n_fiber_per_chan')
        API = dc.find_parameter(params, 'API')

        # Run ratefunc on kwargs and get firing rates for each input
        rates = dc.run_rates_util(ratefunc, params)

        # Check to see if the elements of rates are ndarrays or lists... if they are not lists, we need to put
        # rates inside a list so it can be processed by the list comprehension below
        if type(rates[0]) is not list:
            rates = [rates]

        # Compute partial derivative matrices for rates for AI and then RP
        pdms_AI = [compute_information_curve(x, fs, delta_theta, n_fiber_per_chan, 'AI') for x in rates]
        pdms_AI = np.array(pdms_AI)

        pdms_RP = [compute_information_curve(x, fs, delta_theta, n_fiber_per_chan, 'RP') for x in rates]
        pdms_RP = np.array(pdms_RP)

        return pdms_AI, pdms_RP

    return inner


def compute_information_curve(x, fs, delta_theta, n_fiber_per_chan, _type):
    """ Computes information curves

    Args:
        x (list): list of ndarrays containing firing-rate simulations in shape (n_channel x n_sample). The first
            array should be a firing-rate simulation for baseline parameter values. The following arrays should
            be firing-rate simulations where a single parameter has been incremented by a small amount.
        fs (int): sampling rate in Hz
        delta_theta (ndarray): 1d ndarray containing the increment size for each element of x after the first
        n_fiber_per_chan (array): array containing integers of len n_cf, each element indicates how many fibers
            are theoretically represented by the single corresponding channel in x
        _type (str): either 'AI' or 'RP' for all-information or rate-place

    """
    # Calculate n_param
    n_param = len(x) - 1
    if n_param < 1:
        raise ValueError('There is only one simulation per condition --- ideal observer needs n_param + 1 '
                         'simulations!')
    # Transform from list to ndarray
    x = np.array(x)
    x = np.transpose(x, [1, 0, 2])  # shape: n_cf x (n_param + 1) x n_sample
    # Add small baseline firing rate to avoid issues with zeros and NaNs
    x += 1
    # Construct one ndarray of baseline values and another of incremented values
    baseline = np.tile(x[:, 0, :], [n_param, 1, 1])
    baseline = np.transpose(baseline, [1, 0, 2])  # shape: n_cf x n_param x n_sample
    incremented = x[:, 1:, :]  # shape: n_cf x n_param x n_sample
    if _type == 'AI':
        # Estimate derivative with respect to each parameter
        deriv_estimate = np.transpose(np.transpose((incremented - baseline), [0, 2, 1]) / delta_theta,
                                      [0, 2, 1])  # shape: n_CF x n_param x n_time
        # Normalize the derivatives by the square root of rate
        deriv_norm = np.sqrt(1 / baseline) * deriv_estimate  # shape: n_CF x n_param x n_time
        deriv_matrix = 1 / fs * np.matmul(deriv_norm,
                                          np.transpose(deriv_norm, [0, 2, 1]))  # shape: n_CF x n_param x n_param
        return deriv_matrix
    elif _type == 'RP':
        # Calculate the duration of the response
        t_max = baseline.shape[2] * 1 / fs
        # Average results across time
        baseline = np.mean(baseline, axis=2)
        incremented = np.mean(incremented, axis=2)
        # Estimate derivative with respect to each parameter
        deriv_estimate = (incremented - baseline) / delta_theta
        # Normalize the derivatives by the square root of rate
        deriv_norm = np.sqrt(1 / baseline) * deriv_estimate  # shape: n_CF x n_param
        # Compute derivative matrix
        deriv_norm = np.stack((deriv_norm, deriv_norm), axis=2)
        deriv_matrix = np.matmul(deriv_norm, np.transpose(deriv_norm, [0, 2, 1]))  # shape: n_CF x n_param x n_param
        return deriv_matrix


def simulate_supfigure3_information_curves_fdls(freq, level, model=anf.AuditoryNerveZilany2014, 
                                                model_name='Zilany2014', fs=200e3, 
                                                fs_synapse=100e3, n_cf=40):
    # Define stimulus parameters
    dur = 0.10  # seconds
    dur_ramp = 0.01  # seconds

    # Define model parameters
    cf_low = 0.5*freq
    cf_high = 1.5*freq
    n_cf = n_cf 
    n_fiber_per_chan = round(((np.log10(1.5/0.5)/3)*18000)/n_cf)  # assume ~18k HSR fibers from 0.2 to 20 kHz

    # Encode parameters
    params = si.Parameters(dur=dur, dur_ramp=dur_ramp, fs=fs, n_cf=n_cf, delta_theta=[0.001], API=np.zeros(1),
                           n_fiber_per_chan=n_fiber_per_chan, model_name=model_name, freq=freq, level=level,
                           cf_low=cf_low, cf_high=cf_high, fs_synapse=fs_synapse)

    # Adjust levels to be in dB re: threshold
    params.flatten()
    for ele in params:
        ele['nominal_level'] = ele['level']                                 # encode nominal level (dB re: threshold)
        ele['level'] = adjust_level(ele['freq'], ele['level'], model_name)  # encode actual level (dB SPL)

    # Encode increments
    params.increment({'freq': 0.001})  # increment freq

    # Synthesize stimuli
    synth = sy.PureTone()
    stimuli = synth.synthesize_sequence(params)
    params.add_inputs(stimuli)

    # Construct simulation and run
    sim = model()
    results = sim.run(params, runfunc=get_information_curves(sim.simulate))

    return results


def simulate_supfigure3_information_curves_f0dls(F0, level, model=anf.AuditoryNerveZilany2014, 
                                                 model_name='Zilany2014', fs=200e3, 
                                                 fs_synapse=100e3, n_cf=40):
    # Define stimulus parameters
    dur = 0.10  # seconds
    dur_ramp = 0.01  # seconds

    # Define model parameters
    cf_low = 5*F0
    cf_high = 11*F0
    n_cf = n_cf 
    n_fiber_per_chan = round(((np.log10(11/5)/3)*18000)/n_cf)  # assume ~18k HSR fibers from 0.2 to 20 kHz

    # Encode parameters
    params = si.Parameters(dur=dur, dur_ramp=dur_ramp, fs=fs, n_cf=n_cf, delta_theta=[0.001], API=np.zeros(1),
                           n_fiber_per_chan=n_fiber_per_chan, model_name=model_name, F0=F0, level=level,
                           cf_low=cf_low, cf_high=cf_high, fs_synapse=fs_synapse)

    # Adjust levels to be in dB re: threshold
    params.flatten()
    for ele in params:
        ele['nominal_level'] = ele['level']                                 # encode nominal level (dB re: threshold)
        ele['level'] = adjust_level(ele['F0']*8, ele['level'], model_name)  # encode actual level (dB SPL)

    # Encode increments
    params.increment({'F0': 0.001})  # increment F0

    # Synthesize stimuli
    synth = ISOToneGuest2021_exp1a()
    stimuli = synth.synthesize_sequence(params)
    params.add_inputs(stimuli)

    # Construct simulation and run
    sim = model()
    results = sim.run(params, runfunc=get_information_curves(sim.simulate))

    return results


def simulate_excitation_pattern_f0dls(F0, F0_base, level, model=anf.AuditoryNerveZilany2014, 
                                      model_name='Zilany2014', fs=200e3, 
                                      fs_synapse=100e3, n_cf=40):
    # Define stimulus parameters
    dur = 0.10  # seconds
    dur_ramp = 0.01  # seconds

    # Define model parameters
    cf_low = 5*F0_base
    cf_high = 11*F0_base
    n_cf = n_cf 
    n_fiber_per_chan = round(((np.log10(11/5)/3)*18000)/n_cf)  # assume ~18k HSR fibers from 0.2 to 20 kHz

    # Encode parameters
    params = si.Parameters(dur=dur, dur_ramp=dur_ramp, fs=fs, n_cf=n_cf, delta_theta=[0.001], API=np.zeros(1),
                           n_fiber_per_chan=n_fiber_per_chan, model_name=model_name, F0=F0, level=level,
                           cf_low=cf_low, cf_high=cf_high, fs_synapse=fs_synapse)

    # Adjust levels to be in dB re: threshold
    params.flatten()
    for ele in params:
        ele['nominal_level'] = ele['level']                                 # encode nominal level (dB re: threshold)
        ele['level'] = adjust_level(ele['F0']*8, ele['level'], model_name)  # encode actual level (dB SPL)

    # Synthesize stimuli
    synth = ISOToneGuest2021_exp1a()
    stimuli = synth.synthesize_sequence(params)
    params.add_inputs(stimuli)

    # Construct simulation and run
    sim = model()
    results = sim.run(params)

    return np.mean(results[0], axis=1)


def simulate_excitation_pattern_fdls(freq, freq_base, level, model=anf.AuditoryNerveZilany2014, 
                                      model_name='Zilany2014', fs=200e3, 
                                      fs_synapse=100e3, n_cf=40):
    # Define stimulus parameters
    dur = 0.10  # seconds
    dur_ramp = 0.01  # seconds

    # Define model parameters
    cf_low = 0.5*freq_base
    cf_high = 1.5*freq_base
    n_cf = n_cf 
    n_fiber_per_chan = round(((np.log10(1.5/0.5)/3)*18000)/n_cf)  # assume ~18k HSR fibers from 0.2 to 20 kHz

    # Encode parameters
    params = si.Parameters(dur=dur, dur_ramp=dur_ramp, fs=fs, n_cf=n_cf, delta_theta=[0.001], API=np.zeros(1),
                           n_fiber_per_chan=n_fiber_per_chan, model_name=model_name, freq=freq, level=level,
                           cf_low=cf_low, cf_high=cf_high, fs_synapse=fs_synapse)

    # Adjust levels to be in dB re: threshold
    params.flatten()
    for ele in params:
        ele['nominal_level'] = ele['level']                                 # encode nominal level (dB re: threshold)
        ele['level'] = adjust_level(ele['freq'], ele['level'], model_name)  # encode actual level (dB SPL)

    # Synthesize stimuli
    synth = sy.PureTone()
    stimuli = synth.synthesize_sequence(params)
    params.add_inputs(stimuli)

    # Construct simulation and run
    sim = model()
    results = sim.run(params)
    return np.mean(results[0], axis=1)


def preprocess_information_curves(curves):
    return np.squeeze(curves[0]), np.squeeze(curves[1])


def plot_information_curves(ax, curve, F0):
    cfs = 10**np.linspace(np.log10(F0*5), np.log10(F0*11), n_cf)
    ax.plot(cfs/F0, curve)
    #for ii in range(1, 20):
    #    ax.plot([ii, ii], [0, 10000], linestyle='dashed', color='gray')
    for ii in [6, 7, 8, 9, 10]:
        ax.plot([ii, ii], [0, 10000], linestyle='dashed', color='black')
    ax.set_xlim((4, 12))
    ax.set_xticks([5, 6, 7, 8, 9, 10, 11])
    ax.set_xticklabels([5, 6, 7, 8, 9, 10, 11])
    ax.set_xscale('log')
    ax.set_yscale('log')


def plot_information_curves_fdls(ax, curve, freq):
    cfs = 10**np.linspace(np.log10(freq*0.5), np.log10(freq*1.5), n_cf)
    ax.plot(cfs/freq, curve)
    ax.plot([1, 1], [0, 10000], linestyle='dashed', color='black')
    ax.set_xlim((0.25, 1.75))
    #ax.set_xticks([5, 6, 7, 8, 9, 10, 11])
    #ax.set_xticklabels([5, 6, 7, 8, 9, 10, 11])
    #ax.set_xscale('log')
    ax.set_yscale('log')



def plot_excitation_patterns(ax, curve, F0, color='black', linestyle='solid', label=''):
    cfs = 10**np.linspace(np.log10(F0*5), np.log10(F0*11), n_cf)
    ax.plot(cfs/F0, curve, color=color, linestyle=linestyle, label=label)
    for ii in [6, 7, 8, 9, 10]:
        ax.plot([ii, ii], [0, 10000], linestyle='dashed', color='black', label='_')
    ax.set_xlim((4, 12))
    ax.set_xscale('log')


def plot_excitation_patterns_fdls(ax, curve, freq, color='black', linestyle='solid', label=''):
    cfs = 10**np.linspace(np.log10(freq*0.5), np.log10(freq*1.5), n_cf)
    ax.plot(cfs/freq, curve, color=color, linestyle=linestyle, label=label)
    ax.plot([1, 1], [0, 10000], linestyle='dashed', color='black', label='_')
    ax.set_xlim((0.25, 1.75))
    #ax.set_xticks([5, 6, 7, 8, 9, 10, 11])
    #ax.set_xticklabels([5, 6, 7, 8, 9, 10, 11])
    #ax.set_xscale('log')


def plot_excitation_pattern_difference(ax, curve1, curve2, F0, color):
    cfs = 10**np.linspace(np.log10(F0*5), np.log10(F0*11), n_cf)
    plot_excitation_patterns(ax, curve1, F0, color=color, linestyle='solid')
    plot_excitation_patterns(ax, curve2, F0, color=color, linestyle='dashed', label='_')
    ax.fill_between(cfs/F0, curve1, curve2, curve1 < curve2, color='green')
    ax.fill_between(cfs/F0, curve1, curve2, curve1 > curve2, color='red')


def plot_excitation_pattern_difference_fdls(ax, curve1, curve2, freq, color):
    cfs = 10**np.linspace(np.log10(freq*0.5), np.log10(freq*1.5), n_cf)
    plot_excitation_patterns_fdls(ax, curve1, freq, color=color, linestyle='solid')
    plot_excitation_patterns_fdls(ax, curve2, freq, color=color, linestyle='dashed', label='_')
    ax.fill_between(cfs/freq, curve1, curve2, curve1 < curve2, color='green')
    ax.fill_between(cfs/freq, curve1, curve2, curve1 > curve2, color='red')


# # F0DLs plot
# n_cf = 250
# results_low = [preprocess_information_curves(simulate_supfigure3_information_curves_f0dls(280, level, n_cf=n_cf)[0]) for level in [20, 80]]
# results_high = [preprocess_information_curves(simulate_supfigure3_information_curves_f0dls(1400, level, n_cf=n_cf)[0]) for level in [20, 80]]
# fig, axs = plt.subplots(2, 2, sharey=True, sharex=True)
# for sim in [0, 1]:
#     plot_information_curves(axs[0, 0], results_low[sim][0], 280)
#     plot_information_curves(axs[1, 0], results_low[sim][1], 280)
#     plot_information_curves(axs[0, 1], results_high[sim][0], 1400)
#     plot_information_curves(axs[1, 1], results_high[sim][1], 1400)


# # FDLs plot
# n_cf = 250
# results_low = [preprocess_information_curves(simulate_supfigure3_information_curves_fdls(280*8, level, n_cf=n_cf)[0]) for level in [20, 80]]
# results_high = [preprocess_information_curves(simulate_supfigure3_information_curves_fdls(1400*8, level, n_cf=n_cf)[0]) for level in [20, 80]]
# fig, axs = plt.subplots(2, 2, sharey=True, sharex=True)
# for sim in [0, 1]:
#     plot_information_curves_fdls(axs[0, 0], results_low[sim][0], 280*8)
#     plot_information_curves_fdls(axs[1, 0], results_low[sim][1], 280*8)
#     plot_information_curves_fdls(axs[0, 1], results_high[sim][0], 1400*8)
#     plot_information_curves_fdls(axs[1, 1], results_high[sim][1], 1400*8)

# FDL vs F0DLs information plot
n_cf = 250
results_fdls = [preprocess_information_curves(simulate_supfigure3_information_curves_fdls(1400*8, level, n_cf=n_cf)[0]) for level in [20, 80]]
results_f0dls = [preprocess_information_curves(simulate_supfigure3_information_curves_f0dls(1400, level, n_cf=n_cf)[0]) for level in [20, 80]]
fig, axs = plt.subplots(2, 1, sharey=True, figsize=(4, 5))
for sim, color in zip([0, 1], ['black', 'blue']):
    plot_information_curves_fdls(axs[0], results_fdls[sim][1], 1400)
    plot_information_curves(axs[1], results_f0dls[sim][1], 1400)
axs[0].set_xlabel('CF/frequency')
axs[1].set_xticks([6, 7, 8, 9, 10])
axs[1].set_xticklabels([6, 7, 8, 9, 10])
axs[1].set_xlabel('CF/F0')
axs[0].set_ylabel('Information')
axs[1].set_ylabel('Information')
axs[0].set_ylim((1e-6, 1e0))
axs[1].set_ylim((1e-6, 1e0))
plt.tight_layout()
plt.savefig(os.path.join('plots', 'supfigure_levels_information_patterns.png'), dpi=300)

# FDL vs F0DLs excitation_pattern plot
n_cf = 250
results_fdls = [simulate_excitation_pattern_fdls(1400*8, 1400*8, level, n_cf=n_cf) for level in [20, 80]]
results_fdls_inc = [simulate_excitation_pattern_fdls(1400*2**(1/12)*8, 1400*8, level, n_cf=n_cf) for level in [20, 80]]
results_f0dls = [simulate_excitation_pattern_f0dls(1400, 1400, level, n_cf=n_cf) for level in [20, 80]]
results_f0dls_inc = [simulate_excitation_pattern_f0dls(1400*2**(1/12), 1400, level, n_cf=n_cf) for level in [20, 80]]
fig, axs = plt.subplots(2, 1, sharey=True, figsize=(4, 5))
for sim, color in zip([0, 1], ['blue', 'orange']):
    plot_excitation_pattern_difference_fdls(axs[0], results_fdls[sim], results_fdls_inc[sim], 1400, color=color)
    plot_excitation_pattern_difference(axs[1], results_f0dls[sim], results_f0dls_inc[sim], 1400, color=color)
axs[0].set_ylim((0, 350))
axs[1].set_ylim((0, 350))
axs[0].set_xlabel('CF/frequency')
axs[1].set_xticks([6, 7, 8, 9, 10])
axs[1].set_xticklabels([6, 7, 8, 9, 10])
axs[1].set_xlabel('CF/F0')
axs[0].set_ylabel('Firing rate (sp/s)')
axs[1].set_ylabel('Firing rate (sp/s)')
plt.tight_layout()
plt.savefig(os.path.join('plots', 'supfigure_levels_excitation_patterns.png'), dpi=300)

# FDL vs F0DLs excitation_pattern plot
n_cf = 250
fig, ax = plt.subplots(1, 1, sharey=True, figsize=(4, 2))
for sim, color in zip([0, 1], ['blue', 'orange']):
    plot_excitation_pattern_difference(ax, results_f0dls[sim], results_f0dls_inc[sim], 1400, color=color)
ax.set_ylim((270, 295))
ax.set_xticks([6, 7, 8, 9, 10])
ax.set_xticklabels([6, 7, 8, 9, 10])
ax.set_xlabel('CF/F0')
ax.set_ylabel('Firing rate (sp/s)')
plt.tight_layout()
plt.savefig(os.path.join('plots', 'supfigure_levels_excitation_patterns_zoom.png'), dpi=300)