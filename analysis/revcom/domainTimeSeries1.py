#!/usr/bin/env python3

from science.utility import makeSuperDict
from science.parsing import loadCol

import matplotlib
import matplotlib.pyplot as plt

sims = ['6ZGD_7', '6ZGD_4', '4HFI_7', '4HFI_4']
fancysim = ["Closed, pH 7.0", "Closed, pH 4.0", "Open, pH 7.0", "Open, pH 4.0"]
reps    = [1, 2, 3, 4]
chains  = ['A', 'B', 'C', 'D', 'E']
metrics = ['ecd_twist', 'ecd_spread', 'ecd_upper_spread', 'beta_expansion', 'm2_m1_dist', 'm2_radius', 'm1_kink', 'm1_kink_alt', 'nine_prime_dist', 'nine_prime_pore', 'minus_two_prime_dist', 'minus_two_prime_pore', 'c_loop']
fancy   = ['ECD twist', 'ECD spread', 'Upper ECD Spread', 'beta expansion', 'M1-M2 distance', 'M2 radius', 'M1 kink', 'M1 kink alter', '9\' distance', '9\' pore', '-2\' distance', '-2\' pore', 'loop C']

# Load existing data from ../bloom.

superData = makeSuperDict([sims, reps, chains, metrics, []])

for sim in sims:
    for rep in reps:
        for chain in chains:
            fname = '../bloom/{}_{}_{}.txt'.format(sim, rep, chain)
            for idx in range(0, len(metrics)):
                superData[sim][rep][chain][metrics[idx]] = loadCol(fname, idx + 1, header=0)

# Do plotting.

# Set global font size for figures.
matplotlib.rcParams.update({'font.size': 22})

nrows = 4
ncols = 4

fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(17, 11), dpi=150)

for row in range(0, nrows):
    for col in range(0, ncols):

        subplt = axs[row, col]

        for chain in chains:
            x = superData[sims[col]][reps[row]][chain]["ecd_twist"]
            t = range(0, len(x))

            # t = loadCol(f"{sims[col]}_{reps[row]}_{chain}.txt", col=1)
            # x = loadCol(f"{sims[col]}_{reps[row]}_{chain}.txt", col=2)
            subplt.plot(t, x, linewidth=0.5)

        # Set x-lim and y-lim.
        subplt.set_ylim([-30, -5]) # ECD Twist.
        subplt.set_xlim([0, 1000])

        subplt.text(15, -29, str(row+1))

        # If we're not in the last row, do not show the xticks.
        if row != nrows - 1:
            subplt.set_xticks([])
        else:
            subplt.set_xlabel("Time (ns)")
            subplt.set_xticks([0, 250, 500, 750, 1000])

        # If we're not in the first column, do not show the yticks.
        if col != 0:
            subplt.set_yticks([])
        else:
            subplt.set_ylabel(r'ECD Twist ($\AA$)')
            subplt.set_yticks([-30, -25, -20, -15, -10, -5])

        # Add title to top row.
        if row == 0:
            subplt.set_title(fancysim[col])

fig.tight_layout(pad=0.3)
fig.savefig('ecdtwist.png')
fig.clf()
