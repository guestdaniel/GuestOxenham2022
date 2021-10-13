"""
This script implements the simulations described in Figure 6 of Guest and Oxenham (2021).
"""
import apcmodels.simulation as si
import apcmodels.anf as anf
import apcmodels.decode as dc
from apcmodels.util import save_to_csv
import numpy as np
import os, sys
sys.path.append(os.getcwd())
from util.functions import ISOToneGuest2021, GEOMToneGuest2021, adjust_level
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


def simulate_supfigure3_information_curves_iso(F0, level, model, model_name, fs, n_cf=40):
    """
    Estimates F0 difference limens (FDLs) using ideal observer analysis for a given auditory nerve model. Saves
    the results to disk. This specific harmonic complex tone stimulus used in this simulation is from Guest and
    Oxenham (2021), although no acoustic noise is included in the stimulus.

    Args:
        model: model object from apcmodels.anf
        model_name (str): name of the model
        fs (int): sampling rate in Hz
    """
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
                           cf_low=cf_low, cf_high=cf_high)

    # Adjust levels to be in dB re: threshold
    params.flatten()
    for ele in params:
        ele['nominal_level'] = ele['level']                                 # encode nominal level (dB re: threshold)
        ele['level'] = adjust_level(ele['F0']*8, ele['level'], model_name)  # encode actual level (dB SPL)

    # Encode increments
    params.increment({'F0': 0.001})  # increment F0

    # Synthesize stimuli
    synth = ISOToneGuest2021()
    stimuli = synth.synthesize_sequence(params)
    params.add_inputs(stimuli)

    # Construct simulation and run
    sim = model()
    results = sim.run(params, runfunc=get_information_curves(sim.simulate))

    return results


def simulate_supfigure3_information_curves_geom(F0, level, model, model_name, fs, n_cf=40):
    """
    Estimates F0 difference limens (FDLs) using ideal observer analysis for a given auditory nerve model. Saves
    the results to disk. This specific harmonic complex tone stimulus used in this simulation is from Guest and
    Oxenham (2021), although no acoustic noise is included in the stimulus.

    Args:
        model: model object from apcmodels.anf
        model_name (str): name of the model
        fs (int): sampling rate in Hz
    """
    # Define stimulus parameters
    dur = 0.10  # seconds
    dur_ramp = 0.01  # seconds

    # Define model parameters
    cf_low = 5*F0
    cf_high = 11*F0
    n_cf = n_cf 
    n_fiber_per_chan = round(((np.log10(11/5)/3)*18000)/n_cf)  # assume ~18k HSR fibers from 0.2 to 20 kHz

    # Encode parameters
    params = si.Parameters(dur=dur, dur_ramp=dur_ramp, fs=fs, n_cf=n_cf, delta_theta=[0.001, 0.001], API=np.zeros((2,2)),
                           n_fiber_per_chan=n_fiber_per_chan, model_name=model_name, 
                           F0=lambda: F0*np.random.uniform(0.9, 1.1, 1)[0], 
                           F0_masker=lambda: F0*np.random.uniform(0.9, 1.1, 1)[0],
                           nominal_F0=F0,
                           level=level, cf_low=cf_low, cf_high=cf_high)

    # Adjust levels to be in dB re: threshold
    params.flatten()
    for ele in params:
        ele['nominal_level'] = ele['level']                                 # encode nominal level (dB re: threshold)
        ele['level'] = adjust_level(ele['nominal_F0']*8, ele['level'], model_name)  # encode actual level (dB SPL)

    # Encode increments
    params.increment({'F0': 0.001, 'F0_masker': 0.001})  # increment F0

    # Synthesize stimuli
    synth = GEOMToneGuest2021()
    stimuli = synth.synthesize_sequence(params)
    params.add_inputs(stimuli)

    # Construct simulation and run
    sim = model()
    results = sim.run(params, runfunc=get_information_curves(sim.simulate))

    return results



def preprocess_information_curves(curves):
    return np.squeeze(curves[0]), np.squeeze(curves[1])


def plot_information_curves(ax, curve, F0):
    cfs = 10**np.linspace(np.log10(F0*5), np.log10(F0*11), n_cf)
    ax.plot(cfs, curve)
    for ii in range(1, 20):
        ax.plot([ii*F0, ii*F0], [0, 10000], linestyle='dashed', color='gray')
    ax.set_xlim((4*F0, 12*F0))
    #ax.set_ylim((1e-6, 1e4))
    ax.set_xscale('log')
    ax.set_yscale('log')


def plot_information_curves_fill_between(ax, curve1, curve2, F0):
    # Calculate CFs
    cfs = 10**np.linspace(np.log10(F0*5), np.log10(F0*11), n_cf)
    
    # Plot!
    ax.plot(cfs/F0, curve1, color='black')
    ax.plot(cfs/F0, curve2, color='gray')
    ax.fill_between(cfs/F0, curve1, curve2, color='red', where=curve1 > curve2)
    ax.fill_between(cfs/F0, curve1, curve2, color='green', where=curve1 < curve2)

    # Plot harmonic indicators
    for harm in range(1, 20):
        # Select color
        if harm in [6, 7, 8, 9, 10]:
            color = 'darkgray'
        else:
            color = 'lightgray'
        ax.plot([harm, harm], [0, 10000], 
                linestyle='dashed', color=color)

    # Set limits
    ax.set_xlim((4.5, 11.5))
    ax.set_ylim((1e-6, 1e4))
    ax.set_xscale('log')
    ax.set_yscale('log')


# Get one result
n_cf = 50
iso = preprocess_information_curves(simulate_supfigure3_information_curves_iso(280, 40, anf.AuditoryNerveHeinz2001, 'Heinz2001', fs=500e3, n_cf=n_cf)[0])
geom = preprocess_information_curves(simulate_supfigure3_information_curves_geom(280, 40, anf.AuditoryNerveHeinz2001, 'Heinz2001', fs=500e3, n_cf=n_cf)[0])

# Plot information curves
fig, axs = plt.subplots(2, 1)
# Plot ISO
plot_information_curves(axs[0], iso[0], 280)
plot_information_curves(axs[1], iso[1], 280)
# Plot GEOM
plot_information_curves(axs[0], geom[0][:, 0, 0], 280)
plot_information_curves(axs[1], geom[1][:, 0, 0], 280)

# Plot information curves with fill between
figs, axs = plt.subplots(2, 1, sharex=True)
plot_information_curves_fill_between(axs[0], iso[0], geom[0][:, 0, 0], 280)
plot_information_curves_fill_between(axs[1], iso[1], geom[1][:, 0, 0], 280)

# Note: at this point, we can clearly see that, consistent with threshold predictions, 
# the GEOM stimulus has roughly the same information as the ISO stimulus from the point
# of view of the all-information observer, but substantially *more* information from the 
# point of view of the rate-place observer... Now we should try to visualize an example
# to understand why this is the case!


def generate_excitation_pattern(stimulus, F0=280, F0_base=280, masker_interval=1, 
                                model=anf.AuditoryNerveHeinz2001, model_name='Heinz2001', 
                                level=30, n_cf=40, fs=500e3):
    # Define stimulus parameters
    dur = 0.10  # seconds
    dur_ramp = 0.01  # seconds

    # Define model parameters
    cf_low = 5*F0_base
    cf_high = 11*F0_base
    n_cf = n_cf 

    # Encode parameters
    params = si.Parameters(dur=dur, dur_ramp=dur_ramp, fs=fs, n_cf=n_cf,
                           model_name=model_name, 
                           F0=F0, F0_masker=F0_base*2**(masker_interval/12),
                           nominal_F0=F0,
                           level=level, cf_low=cf_low, cf_high=cf_high)

    # Adjust levels to be in dB re: threshold
    params.flatten()
    for ele in params:
        ele['nominal_level'] = ele['level']                                 # encode nominal level (dB re: threshold)
        ele['level'] = adjust_level(ele['nominal_F0']*8, ele['level'], model_name)  # encode actual level (dB SPL)

    # Synthesize stimuli
    params.add_inputs(stimulus().synthesize_sequence(params))

    # Construct simulation and run
    sim = model()
    results = sim.run(params)

    return np.mean(results[0], axis=1)


def plot_iso_vs_geom_excitation_patterns(F0=280, tar_size=0.3, masker_interval=1, n_cf=200):
    # Set up figure
    fig, axs = plt.subplots(2, 1)

    # Calculate cfs
    cfs = 10**np.linspace(np.log10(F0*5), np.log10(F0*11), n_cf)

    # Plot ISO
    ref = generate_excitation_pattern(ISOToneGuest2021, F0=F0*2**(0/12), n_cf=n_cf)
    inc = generate_excitation_pattern(ISOToneGuest2021, F0=F0*2**(tar_size/12), n_cf=n_cf)
    axs[0].plot(cfs, ref, color='black', linewidth=2)
    axs[0].plot(cfs, inc, color='gray', linewidth=2)
    axs[0].fill_between(cfs, ref, inc, color='red', where=inc < ref)
    axs[0].fill_between(cfs, ref, inc, color='green', where=inc > ref)
    # Plot harmonics
    for harm in [6, 7, 8, 9, 10]:
        axs[0].plot([F0*harm, F0*harm], [0, 150], color='gray', linestyle='dashed')
    # Set limits
    axs[0].set_ylim((0, 150))

    # Calculate negative and positive areas
    x = inc - ref
    x[x < 0] = 0
    pos_area = np.sum(x)
    x = ref - inc
    x[x < 0] = 0
    neg_area = np.sum(x)
    axs[0].set_title("Pos: " + str(pos_area) + " Neg: " + str(neg_area))


    # Plot GEOM
    ref = generate_excitation_pattern(GEOMToneGuest2021, F0=F0*2**(0/12), n_cf=n_cf, masker_interval=masker_interval)
    inc = generate_excitation_pattern(GEOMToneGuest2021, F0=F0*2**(tar_size/12), n_cf=n_cf, masker_interval=masker_interval)
    axs[1].plot(cfs, ref, color='black', linewidth=2)
    axs[1].plot(cfs, inc, color='gray', linewidth=2)
    axs[1].fill_between(cfs, ref, inc, color='red', where=inc < ref)
    axs[1].fill_between(cfs, ref, inc, color='green', where=inc > ref)
    # Plot harmonics
    for harm in [6, 7, 8, 9, 10]:
        axs[1].plot([F0*harm, F0*harm], [0, 150], color='gray', linestyle='dashed')
    # Set limits
    axs[1].set_ylim((0, 150))

    # Calculate negative and positive areas
    x = inc - ref
    x[x < 0] = 0
    pos_area = np.sum(x)
    x = ref - inc
    x[x < 0] = 0
    neg_area = np.sum(x)
    axs[1].set_title("Pos: " + str(pos_area) + " Neg: " + str(neg_area))


plot_iso_vs_geom_excitation_patterns(tar_size=0.3, n_cf=40)