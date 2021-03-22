"""
This script implements the simulations described in Figure TODO of Guest and Oxenham (2021).
"""
import apcmodels.simulation as si
import apcmodels.synthesis as sy
import apcmodels.signal as sg
import apcmodels.anf as anf
from apcmodels.decode import *
import numpy as np
import os, sys
sys.path.append(os.getcwd())
from util.functions import adjust_level
import warnings
warnings.filterwarnings('ignore')
from matplotlib import colors
import matplotlib.pyplot as plt


class ISOToneGuest2021_exp1a(sy.Synthesizer):
    """
    Synthesizes the ISO stimulus in Guest and Oxenham (2021).
    """
    def __init__(self):
        super().__init__(stimulus_name='ISO Tone')

    def synthesize(self, F0, dur=0.350, dur_ramp=0.02, level=None, phase=None, fs=int(48e3), **kwargs):
        """
        Synthesizes a harmonic complex tone composed of components 6-10 of the F0. This is the same stimulus as used
        in Experiment 1a.

        Arguments:
            F0 (float): F0 of the complex tone in Hz
            level (float, ndarray): level per-component of the complex tone in dB SPL, can be either a float (in this
                case, the same level is used for all components) or an ndarray indicating the level of each component
            phase (float, ndarray): phase offset applied to each component of the complex tone in degrees, can be
                either a float (in which case the same phase offset is used for all components) or an ndarray indicating
                the phase offset of each component.
            dur (float): duration in seconds
            dur_ramp (float): duration of raised-cosine ramp in seconds
            fs (int): sampling rate in Hz

        Returns:
            output (array): complex tone stimulus
        """
        # Create array of frequencies, levels, and phases
        freqs = F0*np.arange(1, 21)
        if level is None:
            level = 40*np.ones(len(freqs))  # default to 40 dB SPL per component
        elif isinstance(level, float) or isinstance(level, np.float64) or isinstance(level, int) or isinstance(level, np.int64):
            level = level*np.ones(len(freqs))  # default to 40 dB SPL per component
        if phase is None:
            phase = np.zeros(len(freqs))  # default to sine phase
        elif isinstance(phase, float) or isinstance(phase, np.float64) or isinstance(phase, int) or isinstance(phase, np.int64):
            phase = phase + np.zeros(len(freqs))
        # Synthesize, filter, and ramp complex tone signal
        signal = sg.complex_tone(freqs, level, phase, dur, fs)
        signal = sg.cosine_ramp(signal, dur_ramp, fs)
        # Return
        return signal


def decode_ideal_observer_derivative_mats(ratefunc):
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
            # Compute derivative matrix
            deriv_matrix = 1 / fs * np.matmul(deriv_norm, np.transpose(deriv_norm, [0, 2, 1]))  # shape: n_CF x n_param x n_param
            # Sum across fibers
            deriv_matrix = np.sum(np.transpose(n_fiber_per_chan * np.transpose(deriv_matrix, [1, 2, 0]), [2, 0, 1]), axis=0)
            return deriv_matrix
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
            # Compute derivative matrix
            deriv_norm = np.stack((deriv_norm, deriv_norm), axis=2)
            deriv_matrix = np.matmul(deriv_norm, np.transpose(deriv_norm, [0, 2, 1]))  # shape: n_CF x n_param x n_param
            # Sum across fibers
            deriv_matrix = 0.5 * t_max * np.sum(n_fiber_per_chan * deriv_matrix, axis=0)  # shape: n_param x n_param
            return deriv_matrix

    return inner


def run_rates_util(ratefunc, params):
    """
    Takes inputs and processes each element recursively.

    Arguments:
        ratefunc (function): a function that accepts input and other kwargs and returns model simulations

        params (dict, list): inputs and parameters encoded as a dict or list of dicts. If the input is just a single
            dict, we unpack it and pass it directly to ratefunc. Otherwise, we operate recursively on it.

    Returns:
        output: results of applying ratefunc to each input in params

    """
    # If the input is not a list, just run ratefunc
    output = []
    if type(params) is dict:
        return ratefunc(params)
    # If the input *is* a list, process each input separately
    elif type(params) is list:
        for _input_element in params:
            output.append(run_rates_util(ratefunc, _input_element))
    else:
        raise ValueError('params ought to be a dict or a list')
    return output


def simulate_figure5_f0dls_phase_roving(F0s, levels, model, model_name, fs, n_rep=10):
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
    cf_low = 0.5*F0s
    cf_high = 21*F0s
    n_cf = 80
    n_fiber_per_chan = 40  # TODO: check this value

    # Encode parameters
    delta_theta = [0.001]*21
    API = np.zeros((21, 21))
    API[np.diag_indices(21)] = [0] + [1/360**2]*20
    params = {'dur': dur, 'dur_ramp': dur_ramp, 'fs': fs}  # encode fixed parameters in a dict
    params['phase'] = lambda: np.random.uniform(0, 360, 20)
    params = si.wiggle_parameters(params, 'F0', F0s)  # wiggle frequencies
    params = si.stitch_parameters(params, 'cf_low', cf_low)  # stitch cf_low (each F0 corresponds to a cf_low)
    params = si.stitch_parameters(params, 'cf_high', cf_high)  # stitch cf_high (each F0 corresponds to a cf_high)
    params = si.wiggle_parameters(params, 'level', levels)  # wiggle levels
    params = si.append_parameters(params, ['n_cf', 'delta_theta', 'API', 'n_fiber_per_chan', 'model_name'],
                                  [n_cf, delta_theta, API, n_fiber_per_chan, model_name])  # append other model parameters

    # Adjust levels to be in dB re: threshold
    params = si.flatten_parameters(params)  # flatten out params

    # Encode repeats and increments
    params = si.repeat_parameters(params, n_rep)
    incs = dict()
    incs['F0'] = 0.001
    for ii in range(20):
        key = '(' + str(ii+1) + ')_phase'
        val = np.zeros(20)
        val[ii] = 0.001
        incs[key] = val
    params = si.increment_parameters(params, incs)

    # Synthesize stimuli
    synth = ISOToneGuest2021_exp1a()
    stimuli = synth.synthesize_sequence(params)
    params = si.stitch_parameters(params, '_input', stimuli)

    # Construct simulation and run
    sim = model()
    results = sim.run(params, runfunc=decode_ideal_observer(sim.simulate), parallel=True, hide_progress=False)

    return results


def simulate_figure5_f0dls_level_roving(F0s, levels, model, model_name, fs, n_rep=10):
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
    dur = 0.10  # seconds
    dur_ramp = 0.01  # seconds

    # Define model parameters
    cf_low = 1*F0s
    cf_high = 21*F0s
    n_cf = 80
    n_fiber_per_chan = 40  # TODO: check this value

    # Encode parameters
    delta_theta = [0.001]*21
    API = np.zeros((21, 21))
    API[np.diag_indices(21)] = [0] + [1/6**2] * 20
    params = {'dur': dur, 'dur_ramp': dur_ramp, 'fs': fs}  # encode fixed parameters in a dict
    params = si.wiggle_parameters(params, 'F0', F0s)  # wiggle frequencies
    params = si.stitch_parameters(params, 'cf_low', cf_low)  # stitch cf_low (each F0 corresponds to a cf_low)
    params = si.stitch_parameters(params, 'cf_high', cf_high)  # stitch cf_high (each F0 corresponds to a cf_high)
    params = si.wiggle_parameters(params, 'level', levels)  # wiggle levels
    params = si.append_parameters(params, ['n_cf', 'delta_theta', 'API', 'n_fiber_per_chan', 'model_name'],
                                  [n_cf, delta_theta, API, n_fiber_per_chan, model_name])  # append other model parameters

    # Adjust levels to be in dB re: threshold
    params = si.flatten_parameters(params)  # flatten out params

    # Encode repeats and increments
    params = si.repeat_parameters(params, n_rep)
    incs = dict()
    incs['F0'] = 0.001
    for ii in range(20):
        key = '(' + str(ii+1) + ')_level'
        val = np.zeros(20)
        val[ii] = 0.001
        incs[key] = val
    params = si.increment_parameters(params, incs)

    # Synthesize stimuli
    synth = ISOToneGuest2021_exp1a()
    stimuli = synth.synthesize_sequence(params)
    params = si.stitch_parameters(params, '_input', stimuli)

    # Construct simulation and run
    sim = model()
    results = sim.run(params, runfunc=decode_ideal_observer(sim.simulate), parallel=True, hide_progress=False)
    return results


def plot_derivative_matrices(results_phase, results_level):
    fig, axs = plt.subplots(2, 2)
    for ax, result, idx in zip(np.reshape(axs.T, (4,)), results_phase + results_level, range(4)):
        im = ax.imshow(np.mean(result, axis=0), norm=colors.SymLogNorm(linthresh=0.01, linscale=1, vmin=-1e5, vmax=1e5),
                       cmap='RdBu')
    fig.colorbar(im, ax=np.reshape(axs.T, (4,))[2:4], extend='both')


# Loop through models and calculate FDLs for each model
results_phase = simulate_figure5_f0dls_phase_roving(np.array([300]), np.array([30]), anf.AuditoryNerveHeinz2001, 'Heinz2001', int(200e3))
results_level = simulate_figure5_f0dls_level_roving(np.array([300]), np.array([lambda: np.random.uniform(27, 33, 20)]), anf.AuditoryNerveHeinz2001, 'Heinz2001', int(200e3))

# Plot
plot_derivative_matrices(results_phase[0], results_level[0])

