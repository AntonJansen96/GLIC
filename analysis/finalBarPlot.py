#!/usr/bin/env python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pickle
from science.utility import makeSuperDict

matplotlib.rcParams.update({'font.size': 24})

def stderr(array: list) -> float:
    return np.std(array) / np.sqrt(len(array))

sims = ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']
residues = [13, 14, 26, 31, 32, 35, 49, 55, 67, 69, 75, 82, 86, 88, 91, 97, 104, 115, 122, 127, 136, 145, 147, 153, 154, 161, 163, 177, 178, 181, 185, 222, 235, 243, 272, 277, 282]
data = ['ASP', 'GLU', 'HIS'] + ['E14', 'D31', 'D32', 'E35', 'D49', 'E75'] + ['E26', 'E82', 'D97'] + ['D88', 'D115', 'H127', 'D153', 'D161', 'D178'] + ['D122', 'E177', 'E181'] + ['D185', 'H235', 'E272', 'E282']
# data = ['D13', 'D55', 'E67', 'E69', 'D86', 'D91', 'E104', 'D136', 'D145', 'E147', 'D154', 'E163', 'E222', 'E243', 'H277']

#? LOAD DATA FROM PICKLE OBJECT

protoMean = pickle.load(open('/home/anton/GIT/GLIC/analysis/paperBarPlots.obj', 'rb'))
protoErr = makeSuperDict([sims, residues, 0])

for sim in sims:
    for residue in residues:
        protoErr[sim][residue] = stderr(protoMean[sim][residue])
        protoMean[sim][residue] = np.mean(protoMean[sim][residue])

ncols = 9
nrows = 3

fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(1.8 * ncols, 4.5 * nrows), dpi=300)

row = 0
col = 0

for idx in range(0, len(data)):

    subplt = axs[row, col]

    width = 0.3
    x = np.arange(len([1]))

    if data[idx] in ['ASP', 'GLU', 'HIS']:
        mean1 = [{'ASP': 0.00, 'GLU': 0.00, 'HIS': 0.25}[data[idx]]]
        mean2 = [{'ASP': 0.31, 'GLU': 0.64, 'HIS': 1.00}[data[idx]]]

        subplt.bar(x - 0.5 * width, mean1, width, color='#8856a7')
        subplt.text(x - 0.5 * width - 0.14, mean1[0] + 0.01, f"{mean1[0]:.2f}", size=19)

        subplt.bar(x + 0.5 * width, mean2, width, color='#9ebcda')
        subplt.text(x + 0.5 * width - 0.14, mean2[0] + 0.01, f"{mean2[0]:.2f}", size=19)

    else:
        mean4 = [protoMean['6ZGD_7'][float(data[idx][1:])]]
        serr4 = [protoErr['6ZGD_7'][float(data[idx][1:])]]
        subplt.bar(     x - width * 1.5, mean4, width, color='#8856a7')
        subplt.errorbar(x - width * 1.5, mean4, serr4, color='#8856a7', fmt='none', capsize=7, linewidth=3, markeredgewidth=2)

        mean3 = [protoMean['6ZGD_4'][float(data[idx][1:])]]
        serr3 = [protoErr['6ZGD_4'][float(data[idx][1:])]]
        subplt.bar(     x - width / 2.0, mean3, width, color='#9ebcda')
        subplt.errorbar(x - width / 2.0, mean3, serr3, color='#9ebcda', fmt='none', capsize=7, linewidth=3, markeredgewidth=2)

        mean2 = [protoMean['4HFI_7'][float(data[idx][1:])]]
        serr2 = [protoErr['4HFI_7'][float(data[idx][1:])]]
        subplt.bar(     x + width / 2.0, mean2, width, color='#8856a7', edgecolor='w', lw=0, hatch='//', zorder=13)
        subplt.errorbar(x + width / 2.0, mean2, serr2, color='#8856a7', fmt='none', capsize=7, linewidth=3, markeredgewidth=2, zorder=10)

        mean1 = [protoMean['4HFI_4'][float(data[idx][1:])]]
        serr1 = [protoErr['4HFI_4'][float(data[idx][1:])]]
        subplt.bar(     x + width * 1.5, mean1, width, color='#9ebcda', edgecolor='w', lw=0, hatch='\\\\', zorder=12)
        subplt.errorbar(x + width * 1.5, mean1, serr1, color='#9ebcda', fmt='none', capsize=7, linewidth=3, markeredgewidth=2, zorder=11)

    subplt.set_xticks([])           # Disable x-ticks and set x-label.
    subplt.set_xlabel(data[idx])
    subplt.set_ylim([0, 1.0])       # Set y-lim.

    # Remove boxes around plots.
    for loc in ['top', 'right', 'bottom', 'left']:
        subplt.spines[loc].set_visible(False)

    # If we're not in the first column, do not show the yticks.
    if col != 0:
        subplt.set_yticks([])
    else:
        subplt.set_ylabel('Protonation fraction', size=27)
        subplt.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])

    # Increment row and column indices.
    col += 1
    if col == ncols:
        row += 1
        col = 0

#* REMOVE EMPTY BOXES IF LAST ROW WASN'T FILLED.
if (nrows > 1) and (len(data) % ncols != 0):
    left = ncols - len(data) % ncols

    for ii in range(0, left):
        axs[nrows - 1, ncols - 1 - ii].axis('off')

fig.tight_layout(pad=0.8)
fig.savefig('finalbarplot.png')
fig.clf()
