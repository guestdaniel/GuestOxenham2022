"""
This script implements the simulations described in Figure TODO of Guest and Oxenham (2021).
"""
import apcmodels.simulation as si
import apcmodels.anf as anf
from apcmodels.decode import *
import numpy as np
import os, sys
sys.path.append(os.getcwd())
from util.functions import adjust_level, ISOToneGuest2021_exp1a
import warnings
warnings.filterwarnings('ignore')
from matplotlib import colors
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('QT5Agg')
from matplotlib import patches


def get_io_partial_deriv_gram(ratefunc):
    """
    Returns the partial derivative matrices (i.e., samples of Fisher information matrix) for calculating an ideal
    observer. Implemented as a wrapper that can be applied to any ratefunc that accepts kwargs and returns a single
    firing rate simulation in the standard apcmodels style. Extracted and modified from apcmodels source code.

    Arguments:
        ratefunc (function): a function that accepts **kwargs and returns firing rates for a neural simulation
    """
    def inner(params):
        """
        Runs ratefunc on each input encoded in params, then estimates thresholds based on an ideal observer for a
        particular parameter. This requires some additional information to be encoded in params in the form of an
        a priori information matrix (API) and

        Arguments:
            params (dict): parameters and inputs encoded in a dict or in a list of dicts. % TODO: explain possibilities

        Returns:
            thresholds (ndarray): predicted all-information and rate-place thresholds
        """
        # Pull parameters from encoded list/dict of parameters
        fs = find_parameter(params, 'fs')
        delta_theta = find_parameter(params, 'delta_theta')
        n_fiber_per_chan = find_parameter(params, 'n_fiber_per_chan')
        API = find_parameter(params, 'API')

        # Run ratefunc on kwargs and get firing rates for each input
        rates = run_rates_util(ratefunc, params)

        # Check to see if the elements of rates are ndarrays or lists... if they are not lists, we need to put
        # rates inside a list so it can be processed by the list comprehension below
        if type(rates[0]) is not list:
            rates = [rates]

        # Compute partial derivative matrices for rates for AI and then RP
        pdms_AI = [compute_partial_derivative_matrix(x, fs, delta_theta, n_fiber_per_chan, 'AI') for x in rates]
        pdms_AI = np.array(pdms_AI)

        pdms_RP = [compute_partial_derivative_matrix(x, fs, delta_theta, n_fiber_per_chan, 'RP') for x in rates]
        pdms_RP = np.array(pdms_RP)

        # Return ideal observer results
        return pdms_AI, pdms_RP

    def compute_partial_derivative_matrix(x, fs, delta_theta, n_fiber_per_chan, _type):
        """
        Given one list of simulations, computes a partial derivative matrix as in Siebert (1972).

        Arguments:
            x (list): list of ndarrays containing firing-rate simulations in shape (n_channel x n_sample). The first
                array should be a firing-rate simulation for baseline parameter values. The following arrays should
                be firing-rate simulations where a single parameter has been incremented by a small amount.

            fs (int): sampling rate in Hz

            delta_theta (ndarray): 1d ndarray containing the increment size for each element of x after the first
            n_fiber_per_chan (array): array containing integers of len n_cf, each element indicates how many fibers
                are theoretically represented by the single corresponding channel in x

            _type (str): either 'AI' or 'RP' for all-information or rate-place

        Returns:

        """
        # Calculate n_param
        n_param = len(x)-1
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
            deriv_estimate = np.transpose(np.transpose((incremented - baseline), [0, 2, 1]) / delta_theta, [0, 2, 1])  # shape: n_CF x n_param x n_time
            # Normalize the derivatives by the square root of rate
            deriv_norm = np.sqrt(1 / baseline) * deriv_estimate  # shape: n_CF x n_param x n_time
            # Create empty matrix
            deriv_gram = np.zeros((deriv_norm.shape[0], deriv_norm.shape[1], deriv_norm.shape[1], deriv_norm.shape[2]))
            for ii in range(deriv_norm.shape[1]):
                for jj in range(deriv_norm.shape[1]):
                    if jj == ii:
                        deriv_gram[:, ii, jj, :] = deriv_norm[:, ii, :]
                    else:
                        deriv_gram[:, ii, jj, :] = np.multiply(deriv_norm[:, ii, :], deriv_norm[:, jj, :])
            return deriv_gram
        elif _type == 'RP':
            # Calculate the duration of the response
            t_max = baseline.shape[2] * 1/fs
            # Average results across time
            baseline = np.mean(baseline, axis=2)
            incremented = np.mean(incremented, axis=2)
            # Estimate derivative with respect to each parameter
            deriv_estimate = (incremented - baseline)/delta_theta
            # Normalize the derivatives by the square root of rate
            deriv_norm = np.sqrt(1 / baseline) * deriv_estimate  # shape: n_CF x n_param
            # Create empty matrix
            deriv_gram = np.zeros((deriv_norm.shape[0], deriv_norm.shape[1], deriv_norm.shape[1]))
            for ii in range(deriv_norm.shape[1]):
                for jj in range(deriv_norm.shape[1]):
                    if jj == ii:
                        deriv_gram[:, ii, jj] = deriv_norm[:, ii]
                    else:
                        deriv_gram[:, ii, jj] = np.multiply(deriv_norm[:, ii], deriv_norm[:, jj])
            return deriv_gram

    return inner


def simulate_figure5_f0dls(F0s, levels, model, model_name, fs, extractfunc=get_io_partial_deriv_gram,
                           cf_low=5, cf_high=11, n_cf=40, n_rep=5):
    """
    Estimates F0 difference limens (F0DLs) using ideal observer analysis for a given auditory nerve model. Saves
    the results to disk. The simulations include phase randomization (specifically, each component from 6-10 has its
    starting phase randomized on each trial).

    Arguments:
        model: model object from apcmodels.anf

        model_name (str): name of the model

        fs (int): sampling rate in Hz

        n_rep (int): number of repetitions to perform for each simulation
    """
    # Define stimulus parameters
    dur = 0.10  # seconds
    dur_ramp = 0.01  # seconds

    # Define model parameters
    cf_low = cf_low*F0s
    cf_high = cf_high*F0s
    n_fiber_per_chan = 40  # TODO: check this value

    # Encode parameters
    delta_theta = [0.001]
    API = np.zeros((1, 1))
    params = {'dur': dur, 'dur_ramp': dur_ramp, 'fs': fs}  # encode fixed parameters in a dict
    params = si.wiggle_parameters(params, 'F0', F0s)  # wiggle frequencies
    params = si.stitch_parameters(params, 'cf_low', cf_low)  # stitch cf_low (each F0 corresponds to a cf_low)
    params = si.stitch_parameters(params, 'cf_high', cf_high)  # stitch cf_high (each F0 corresponds to a cf_high)
    params = si.wiggle_parameters(params, 'level', levels)  # wiggle levels
    params = si.append_parameters(params, ['n_cf', 'delta_theta', 'API', 'n_fiber_per_chan', 'model_name'],
                                  [n_cf, delta_theta, API, n_fiber_per_chan, model_name])  # append other model parameters

    # Adjust levels to be in dB re: threshold
    params = si.flatten_parameters(params)  # flatten out params
    for ele in params:
        ele['nominal_level'] = ele['level']  # save the nominal level (dB re: threshold)
        ele['level'] = adjust_level(ele['F0']*8, ele['level'], model_name)  # encode the actual level (dB SPL)

    # Encode repeats and increments
    params = si.repeat_parameters(params, n_rep)
    params = si.increment_parameters(params, {'F0': 0.001})

    # Synthesize stimuli
    synth = ISOToneGuest2021_exp1a()
    stimuli = synth.synthesize_sequence(params)
    params = si.stitch_parameters(params, '_input', stimuli)

    # Construct simulation and run
    sim = model()
    results = sim.run(params, runfunc=extractfunc(sim.simulate))

    return results


def simulate_figure5_f0dls_level_roving(F0s, levels, model, model_name, fs, extractfunc=get_io_partial_deriv_gram, n_rep=5,
                                        cf_low=5, cf_high=11, n_cf=40):
    """
    Estimates F0 difference limens (F0DLs) using ideal observer analysis for a given auditory nerve model. Saves
    the results to disk. Includes simulations of level roving.

    Arguments:
        model: model object from apcmodels.anf

        model_name (str): name of the model

        fs (int): sampling rate in Hz

        n_rep (int): number of repetitions to perform for each simulation
    """
    # Define stimulus parameters
    nominal_levels = [20, 30, 40]
    dur = 0.10  # seconds
    dur_ramp = 0.01  # seconds

    # Define model parameters
    cf_low = cf_low*F0s
    cf_high = cf_high*F0s
    n_fiber_per_chan = 40  # TODO: check this value

    # Encode parameters
    delta_theta = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001]
    API = np.zeros((6, 6))
    API[np.diag_indices(6)] = [0, 1/6**2, 1/6**2, 1/6**2, 1/6**2, 1/6**2]
    params = {'dur': dur, 'dur_ramp': dur_ramp, 'fs': fs}  # encode fixed parameters in a dict
    params = si.wiggle_parameters(params, 'F0', F0s)  # wiggle frequencies
    params = si.stitch_parameters(params, 'cf_low', cf_low)  # stitch cf_low (each F0 corresponds to a cf_low)
    params = si.stitch_parameters(params, 'cf_high', cf_high)  # stitch cf_high (each F0 corresponds to a cf_high)
    params = si.wiggle_parameters_parallel(params, ['level', 'nominal_level'], [levels, nominal_levels])  # wiggle levels
    params = si.append_parameters(params, ['n_cf', 'delta_theta', 'API', 'n_fiber_per_chan', 'model_name'],
                                  [n_cf, delta_theta, API, n_fiber_per_chan, model_name])  # append other model parameters

    # Adjust levels to be in dB re: threshold
    params = si.flatten_parameters(params)  # flatten out params

    # Encode repeats and increments
    params = si.repeat_parameters(params, n_rep)
    params = si.increment_parameters(params, {'F0': 0.001,                                  # increment F0
                                              '(1)_level': np.array([0.001, 0, 0, 0, 0]),   # increment level of H6
                                              '(2)_level': np.array([0, 0.001, 0, 0, 0]),   # increment level of H7
                                              '(3)_level': np.array([0, 0, 0.001, 0, 0]),   # increment level of H8
                                              '(4)_level': np.array([0, 0, 0, 0.001, 0]),   # increment level of H9
                                              '(5)_level': np.array([0, 0, 0, 0, 0.001])})  # increment level of H10

    # Adjust levels to be in dB re: threshold
    for ele in params:   # loop through elements of params
        for x in ele:    # loop through repeats of params
            for y in x:  # loop through increments of params
                y['level'] = adjust_level(y['F0']*8, y['level'], model_name)  # encode the actual level (dB SPL)

    # Synthesize stimuli
    synth = ISOToneGuest2021_exp1a()
    stimuli = synth.synthesize_sequence(params)
    params = si.stitch_parameters(params, '_input', stimuli)

    # Construct simulation and run
    sim = model()
    results = sim.run(params, runfunc=extractfunc(sim.simulate), parallel=True, hide_progress=False)
    return results


# Get and plot partial derivative "grams"
results_level_roving = simulate_figure5_f0dls_level_roving(np.array([280]), np.array([30]), anf.AuditoryNerveHeinz2001, 'Heinz2001', int(200e3), cf_low=4, cf_high=14, n_cf=200)
cfs = 10**np.linspace(np.log10(4*280), np.log10(14*280), 200)

# Define function to plot data
def plot_pretty_stuff(yaxis=True):
    # Set yscale
    plt.yscale('symlog', linthreshy=0.01)
    plt.ylim([-1e0, 1e0])
    plt.xlim(cfs[0], cfs[-1])

    # Plot markers
    plt.plot([5 * 280, 5 * 280], [-1e1, 1e1], 'k')
    plt.plot([11 * 280, 11 * 280], [-1e1, 1e1], 'k')
    for harm in [6, 7, 8, 9, 10]:
        plt.plot([harm * 280, harm * 280], [-1e1, 1e1], color='gray', linestyle='dotted')

    # Plot transparent rectangles
    rect = patches.Rectangle((0, -1e1), 5 * 280, 2e1, linewidth=1, edgecolor='none', facecolor='white', alpha=0.8,
                             zorder=5)
    plt.gca().add_patch(rect)
    rect = patches.Rectangle((11 * 280, -1e1), 5 * 280, 2e1, linewidth=1, edgecolor='none', facecolor='white',
                             alpha=0.8, zorder=5)
    plt.gca().add_patch(rect)

    # Labels
    plt.ylabel('Fisher information')
    plt.xlabel('CF (Hz)')

    # Disable yaxis
    plt.gca().get_yaxis().set_visible(yaxis)

plt.figure(figsize=(8, 3.5))

plt.subplot(1, 2, 1)
plt.plot(cfs, results_level_roving[0][1][0, :, 0, 0]**2, color='black', linewidth=3)
plt.plot(cfs, results_level_roving[0][1][0, :, 0, 1], color='lightcoral')
plt.plot(cfs, results_level_roving[0][1][0, :, 0, 5], color='darkred')
plot_pretty_stuff()
plt.legend(['F0', 'Level x F0 (H6)', 'Level x F0 (H10)'], framealpha=1)

plt.subplot(1, 2, 2)
plt.plot(cfs, results_level_roving[0][1][0, :, 0, 0]**2, color='black', linewidth=3)
plt.plot(cfs, results_level_roving[0][1][0, :, 1, 1]**2, color='lightsteelblue')
plt.plot(cfs, results_level_roving[0][1][0, :, 5, 5]**2, color='darkblue')
plot_pretty_stuff(yaxis=False)
plt.legend(['F0', 'Level (H6)', 'Level (H10)'], framealpha=1)

plt.tight_layout()