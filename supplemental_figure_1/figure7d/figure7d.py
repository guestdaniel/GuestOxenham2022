import numpy as np
import matplotlib.pyplot as plt
import os, sys
sys.path.append(os.getcwd())
import util as cfg
from util.functions import ComplexToneCedolin2005


def plot_isi_histogram(ax, ISIs, F0, first):
    """
    Plots an interspike-interval histograms from a list of interspike intervals recorded in a simulation.

    Arguments:
        ax (axis): axis object on which this figure should be plotted
        ISIs (np.ndarray): array of interspike intervals, of shape (n_interval, )
        F0 (float): F0 at which the interspike intervals were simulated, used to plot vertical lines indicating
            intervals of the F0
        first (bool): bool indicating whether or not this is the first plot in a column of plots, if True then we
            don't label the x-axis
    """
    # Calculate histogram
    hist_low, edges_low = np.histogram(ISIs, bins=2200, range=(0, 0.22))
    # Plot dashed lines at F0 intervals
    for ii in range(1, 50):
        ax.plot([ii * 1 / F0 * 1000, ii * 1 / F0 * 1000], [0, 14000], color='gray', linestyle='dashed',
                linewidth=0.5)
    # Plot histogram
    ax.plot(edges_low[1:] * 1000, hist_low)
    ax.set_xlim((0, 25))
    if first:
        ax.get_xaxis().set_visible(False)
    else:
        ax.set_xlabel('Interval (s)')
    ax.set_ylabel('Interval count')
    ax.set_title(str(F0) + ' Hz')
    ax.set_ylim((0, 14000))

# Create figure and saxes
fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(4, 4))
ISIs = np.load('supplemental_figure_1/figure7d/isi_' + str(320) + '.npy')
plot_isi_histogram(ax[0], ISIs, 320, True)
ISIs = np.load('supplemental_figure_1/figure7d/isi_' + str(880) + '.npy')
plot_isi_histogram(ax[1], ISIs, 880, False)

# Save plot
plt.savefig('plots/fig7d.png')