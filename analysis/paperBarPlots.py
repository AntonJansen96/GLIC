#!/usr/bin/env python3

import pickle
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from science.utility import makeSuperDict
from science.parsing import loadxvg
from science.cphmd import getLambdaFileIndices
from science.cphmd import protonation
from science.cphmd import deprotonation

matplotlib.rcParams.update({'font.size': 20})

letters  = ['D', 'E', 'E', 'D', 'D', 'E', 'D', 'D', 'E', 'E', 'E', 'E', 'D', 'D', 'D', 'D', 'E', 'D', 'D', 'H', 'D', 'D', 'E', 'D', 'D', 'D', 'E', 'E', 'D', 'E', 'D', 'E', 'H', 'E', 'E', 'H', 'E']
residues = [13, 14, 26, 31, 32, 35, 49, 55, 67, 69, 75, 82, 86, 88, 91, 97, 104, 115, 122, 127, 136, 145, 147, 153, 154, 161, 163, 177, 178, 181, 185, 222, 235, 243, 272, 277, 282]

ecd      = [26, 67, 69, 75, 82, 86, 88, 104, 13, 14, 49, 55, 91, 97, 127, 136, 145, 147, 153, 154, 163, 177, 178, 181, 185]
ecdtmd   = [31, 32, 35, 115, 161, 122, 243]
tmd      = [222, 235, 272, 277, 282]

sims     = ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']
reps     = [1, 2, 3, 4]
chains   = ['A', 'B', 'C', 'D', 'E']

if not os.path.isfile('paperBarPlots.obj'):

    protoMean = makeSuperDict([sims, residues, []])

    for residue in residues:
        lambdaIndices = getLambdaFileIndices('../sims/4HFI_4/01/CA.pdb', residue)

        for sim in sims:
            for rep in reps:
                for idx in lambdaIndices:
                    print(residue, sim, rep, f'cphmd-coord-{idx}.xvg')

                    x = loadxvg(f'../sims/{sim}/{rep:02d}/cphmd-coord-{idx}.xvg', dt=1000, b=0)[1]

                    if residue in [127, 235, 277]:
                        protoMean[sim][residue].append(deprotonation(x))
                    else:
                        protoMean[sim][residue].append(protonation(x))

    pickle.dump(protoMean, open('paperBarPlots.obj', 'wb'))

else:
    protoMean = pickle.load(open('paperBarPlots.obj', 'rb'))

protoErr = makeSuperDict([sims, residues, 0])

for sim in sims:
    for residue in residues:
        array = protoMean[sim][residue]
        protoErr[sim][residue]  = np.std(array) / np.sqrt(len(array))
        protoMean[sim][residue] = np.mean(array)

#? We have now processed the data. Now it's time to template the multiplot.

data  = residues
nrows = 3

# Initialize multifigure object.
ncols    = int(np.ceil(len(data) / nrows))
fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20, 11), dpi=100)

row = 0
col = 0
for idx in range(0, len(data)):

    # Select subplt object.
    subplt = axs[row, col]

    # General plotting stuff.
    width = 0.2
    x = np.arange(len([1]))

    # 6ZGD_7
    mean4 = [protoMean['6ZGD_7'][residues[idx]]]
    serr4 = [protoErr['6ZGD_7'][residues[idx]]]
    subplt.bar(     x - width * 1.5, mean4, width, color='#8856a7')
    subplt.errorbar(x - width * 1.5, mean4, serr4, color='#8856a7', fmt='none', capsize=6, linewidth=2)

    # 6ZGD_4
    mean3 = [protoMean['6ZGD_4'][residues[idx]]]
    serr3 = [protoErr['6ZGD_4'][residues[idx]]]
    subplt.bar(     x - width / 2.0, mean3, width, color='#9ebcda')
    subplt.errorbar(x - width / 2.0, mean3, serr3, color='#9ebcda', fmt='none', capsize=6, linewidth=2)

    # 4HFI_7
    mean2 = [protoMean['4HFI_7'][residues[idx]]]
    serr2 = [protoErr['4HFI_7'][residues[idx]]]
    subplt.bar(     x + width / 2.0, mean2, width, color='#8856a7', edgecolor='w', lw=1, hatch='//')
    subplt.errorbar(x + width / 2.0, mean2, serr2, color='#8856a7', fmt='none', capsize=6, linewidth=2)

    # 4HFI_4
    mean1 = [protoMean['4HFI_4'][residues[idx]]]
    serr1 = [protoErr['4HFI_4'][residues[idx]]]
    subplt.bar(     x + width * 1.5, mean1, width, color='#9ebcda', edgecolor='w', lw=1, hatch='\\\\')
    subplt.errorbar(x + width * 1.5, mean1, serr1, color='#9ebcda', fmt='none', capsize=6, linewidth=2)

    # Disable x-ticks and set x-label.
    subplt.set_xticks([])
    subplt.set_xlabel(f'{letters[idx]}{residues[idx]}')

    # Set y-lim.
    subplt.set_ylim([0, 1.1])

    # Remove boxes around plots.
    subplt.spines['top'].set_visible(False)
    subplt.spines['right'].set_visible(False)
    subplt.spines['bottom'].set_visible(False)
    subplt.spines['left'].set_visible(False)

    # If we're not in the first column, do not show the yticks.
    # Else, do show them.
    if col != 0:
        subplt.set_yticks([])
    else:
        subplt.set_ylabel('Protonation')
        subplt.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])

    # This is simply for the rows and columns.
    col += 1
    if col == ncols:
        col = 0
        row += 1

# Remove the last two empty boxes in third row.
axs[2, 11].axis('off')
axs[2, 12].axis('off')

# Proxy artist for legend
subplt.bar(0, 0, 0, color='#8856a7', label='closed, pH 7')
subplt.bar(0, 0, 0, color='#9ebcda', label='closed, pH 4')
subplt.bar(0, 0, 0, color='#8856a7', label='open, pH 7', edgecolor='w', lw=1, hatch='//')
subplt.bar(0, 0, 0, color='#9ebcda', label='open, pH 4', edgecolor='w', lw=1, hatch='\\\\')

fig.legend(loc=[0.85, 0.06], prop={'size': 21})
fig.tight_layout()
fig.savefig('allproto.png')
os.system('convert allproto.png -trim allproto.png')
fig.clf()
