#!/usr/bin/env python3

#? IMPORTS

import pickle
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from science.utility import makeSuperDict
from science.cphmd import theoreticalMicropKa

#? PARAMETERS AND GLOBAL STUFF

matplotlib.rcParams.update({'font.size': 24})

def stderr(array: list) -> float:
    return np.std(array) / np.sqrt(len(array))

sims = ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']
letters  = ['D', 'E', 'E', 'D', 'D', 'E', 'D', 'D', 'E', 'E', 'E', 'E', 'D', 'D', 'D', 'D', 'E', 'D', 'D', 'H', 'D', 'D', 'E', 'D', 'D', 'D', 'E', 'E', 'D', 'E', 'D', 'E', 'H', 'E', 'E', 'H', 'E']
residues = [13, 14, 26, 31, 32, 35, 49, 55, 67, 69, 75, 82, 86, 88, 91, 97, 104, 115, 122, 127, 136, 145, 147, 153, 154, 161, 163, 177, 178, 181, 185, 222, 235, 243, 272, 277, 282]

#? LOAD DATA FROM PICKLE OBJECT

protoMean = pickle.load(open('paperBarPlots.obj', 'rb'))
protoErr = makeSuperDict([sims, residues, 0])

for sim in sims:
    for residue in residues:
        array = protoMean[sim][residue]
        protoErr[sim][residue] = stderr(array)
        protoMean[sim][residue] = np.mean(array)

#? MAKE RESIDUE SELECTION

cat1 = []  # No pH or state dependence (differences less than 1 pKa unit at pH = 4.0).
cat2 = []  # pH dependence (difference more than 1 pKa unit at pH = 4.0).
cat3 = []  # state dependence (differences more than 1 pKa unit between structures at same pH).
cat4 = []  # both state and pH-dependence.

cutoff = 1  # pKa unit(s) deviation

for residue in residues:

    resChar = letters[residues.index(residue)]

    if resChar == 'D':
        macropKa = 3.65
    if resChar == 'E':
        macropKa = 4.25
    if resChar == 'H':
        macropKa = 6.53

    for sim in sims:

        if sim in ['4HFI_4', '6ZGD_4']:
            pH = 4.0
        if sim in ['4HFI_7', '6ZGD_7']:
            pH = 7.0
            continue

        protonation = protoMean[sim][residue]
        micropKa = theoreticalMicropKa(pH, protonation)

        if (micropKa < (macropKa - cutoff)) or (micropKa > (macropKa + cutoff)):
            print(residue, sim, protonation, f'{macropKa:.2f}', f'{micropKa:.2f}')
            cat2.append(residue)
            break

for residue in residues:
    if residue not in cat2:
        cat1.append(residue)

# MAKE SELECTION FOR STATE-DEPENDENCE
for residue in residues:
    cutoff = 1
    A = theoreticalMicropKa(4.0, protoMean['4HFI_4'][residue])
    B = theoreticalMicropKa(7.0, protoMean['4HFI_7'][residue])
    C = theoreticalMicropKa(4.0, protoMean['6ZGD_4'][residue])
    D = theoreticalMicropKa(7.0, protoMean['6ZGD_7'][residue])

    # SELECT ON LOW-PH STRUCTURE DIFFERENCES
    if (A > C + cutoff) or (A < C - cutoff):
        cat3.append(residue)

    # SELECT ALSO ON HIGH-PH STRUCTURE DIFFERENCE
    # elif (B > D + cutoff) or (B < D - cutoff):
    #     if protoMean['4HFI_7'][residue] < 0.01 or protoMean['6ZGD_7'][residue] < 0.01:
    #         continue
    #     cat3.append(residue)

# MAKE SELECTION FOR BOTH PH AND STATE DEPENDENCE
for residue in residues:
    if residue in cat2 and residue in cat3:
        cat4.append(residue)

for residue in cat3:
    if residue in cat2:
        cat2.pop(cat2.index(residue))

# cat3.remove(178)

#! Remove histidines and add them back manually
for resid in [127, 235, 277]:
    for cat in [cat1, cat2, cat3]:
        try:
            cat.remove(resid)
        except ValueError:
            pass
cat1.append(127)
cat2.append(235)
cat2.append(277)

#? MAIN PLOTTING CODE.

for case in ['cat1', 'cat2', 'cat3', 'cat4']:

    if case == 'cat1':
        data  = cat1
        nrows = 1

    if case == 'cat2':
        data  = cat2
        nrows = 2

    if case == 'cat3':
        data  = cat3
        nrows = 1

    if case == 'cat4':
        data  = cat4
        nrows = 1

    outname = case
    data = sorted(data)
    if not data:
        continue

    # Initialize multifigure object.
    ncols    = int(np.ceil(len(data) / nrows))

    if case == 'cat2':
        ncols = 11

    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(1.8 * ncols, 4.5 * nrows), dpi=300)

    row = 0
    col = 0

    # Loop for creating the subfigures.
    for idx in range(0, len(data)):

        # Select subplt object.
        if nrows > 1:
            subplt = axs[row, col]
        else:
            subplt = axs[col]

        # The (recurring) bar plot block.
        width = 0.3
        x = np.arange(len([1]))

        # subplt.set_figheight(15)
        # subplt.set_figwidth(15)

        mean4 = [protoMean['6ZGD_7'][data[idx]]]
        serr4 = [protoErr['6ZGD_7'][data[idx]]]
        subplt.bar(     x - width * 1.5, mean4, width, color='#8856a7')
        subplt.errorbar(x - width * 1.5, mean4, serr4, color='#8856a7', fmt='none', capsize=7, linewidth=3, markeredgewidth=2)
        # subplt.text(    x - width * 1.5 - 0.15, (mean4[0] + serr4[0]) + 0.01, f'{round(mean4[0], 2):.2f}', size=8)

        mean3 = [protoMean['6ZGD_4'][data[idx]]]
        serr3 = [protoErr['6ZGD_4'][data[idx]]]
        subplt.bar(     x - width / 2.0, mean3, width, color='#9ebcda')
        subplt.errorbar(x - width / 2.0, mean3, serr3, color='#9ebcda', fmt='none', capsize=7, linewidth=3, markeredgewidth=2)
        # subplt.text(    x - width / 2.0 - 0.15, (mean3[0] + serr3[0]) + 0.01, f'{round(mean3[0], 2):.2f}', size=8)

        mean2 = [protoMean['4HFI_7'][data[idx]]]
        serr2 = [protoErr['4HFI_7'][data[idx]]]
        subplt.bar(     x + width / 2.0, mean2, width, color='#8856a7', edgecolor='w', lw=0, hatch='//', zorder=13)
        subplt.errorbar(x + width / 2.0, mean2, serr2, color='#8856a7', fmt='none', capsize=7, linewidth=3, markeredgewidth=2, zorder=10)
        # subplt.text(    x + width / 2.0 - 0.15, (mean2[0] + serr2[0]) + 0.01, f'{round(mean2[0], 2):.2f}', size=8)

        mean1 = [protoMean['4HFI_4'][data[idx]]]
        serr1 = [protoErr['4HFI_4'][data[idx]]]
        subplt.bar(     x + width * 1.5, mean1, width, color='#9ebcda', edgecolor='w', lw=0, hatch='\\\\', zorder=12)
        subplt.errorbar(x + width * 1.5, mean1, serr1, color='#9ebcda', fmt='none', capsize=7, linewidth=3, markeredgewidth=2, zorder=11)
        # subplt.text(    x + width * 1.5 - 0.15, (mean1[0] + serr1[0]) + 0.01, f'{round(mean1[0], 2):.2f}', size=8)

        # Disable x-ticks and set x-label.
        subplt.set_xticks([])
        letter = letters[residues.index(data[idx])]
        subplt.set_xlabel(f'{letter}{data[idx]}')

        # Set y-lim.
        subplt.set_ylim([0, 1.0])

        # Remove boxes around plots.
        for loc in ['top', 'right', 'bottom', 'left']:
            subplt.spines[loc].set_visible(False)

        # If we're not in the first column, do not show the yticks.
        if col != 0:
            subplt.set_yticks([])
        else:
            subplt.set_ylabel('Protonation')
            subplt.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])

        # Increment row and column indices.
        col += 1
        if col == ncols:
            row += 1
            col = 0

    #* REMOVE EMPTY BOXES IF LAST ROW WASN'T FILLED.
    if (nrows > 1) and (len(data) % ncols != 0):

        left = ncols - len(data) % ncols
        # print(len(data), ncols, left)  # debug

        for ii in range(0, left):
            axs[nrows - 1, ncols - 1 - ii].axis('off')

    #* LEGEND
    # subplt.bar(0, 0, 0, color='#8856a7', label='closed, pH 7')
    # subplt.bar(0, 0, 0, color='#9ebcda', label='closed, pH 4')
    # subplt.bar(0, 0, 0, color='#8856a7', label='open, pH 7', edgecolor='w', lw=1, hatch='//')
    # subplt.bar(0, 0, 0, color='#9ebcda', label='open, pH 4', edgecolor='w', lw=1, hatch='\\\\')
    # fig.legend(loc=[0.85, 0.06], prop={'size': 20})

    #* SAVE
    fig.tight_layout(pad=0.8)
    fig.savefig(f'{outname}.png')
    fig.savefig(f'{outname}.eps')
    fig.clf()
