import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os, sys
sys.path.append(os.getcwd())
import util as cfg


# Load in histograms and calculate edges
histograms = np.load('supplemental_figure_1/figure7e/isi_histograms.npy')
neural_harm_nums = np.load('supplemental_figure_1/figure7e/neural_harm_nums.npy')
edges = np.linspace(start=0, stop=20, num=2200)

# Create plot
fig, ax = plt.subplots(1, figsize=(8, 3))
plt.pcolormesh(edges, neural_harm_nums, histograms, cmap='gray_r')
plt.xlabel('Lag time (ms)')
plt.ylabel('CF (Hz)')
plt.xlim((0, 8))
plt.xlabel('Periods of Lower F0')
plt.ylabel('Neural Harmonic Number (CF/F0)')
# Highlight peaks at F0 multiples
for ii in range(10):
    if ii==0:
        continue
    rect = patches.Rectangle((ii-0.05, 2.2), 0.10, 3.8, linewidth=0.75, edgecolor='r', linestyle='--', facecolor='none')
    #ax.add_patch(rect)
    rect = patches.Rectangle(((ii*9/11)-0.05, 2.2), 0.10, 3.8, linewidth=0.75, edgecolor='b', linestyle='--', facecolor='none')
    #ax.add_patch(rect)
# Annotate arrows
ax.annotate("1/F0", xy=(1, 5), xytext=(1.2, 5.2), arrowprops=dict(arrowstyle='->', color='red'), color='red')
ax.annotate("1/F0", xy=(1/(11/9), 5), xytext=(1/(11/9) - 0.6, 5.2), arrowprops=dict(arrowstyle='->', color='blue'), color='blue')
ax.annotate("6/F0", xy=(6, 5), xytext=(6.2, 5.2), arrowprops=dict(arrowstyle='->', color='red'), color='red')
ax.annotate("7/F0", xy=(7/(11/9), 5), xytext=(7/(11/9) - 0.6, 5.2), arrowprops=dict(arrowstyle='->', color='blue'), color='blue')
# Save to disk
plt.tight_layout()
plt.savefig('plots/fig7e.png')