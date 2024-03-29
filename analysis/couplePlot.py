#!/usr/bin/env python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os

from science.utility import makeSuperDict
from science.parsing import loadCol

regul = [['E26', 'V79p'], ['E26', 'N80p'], ['E26', 'V81p'], ['E82', 'T36'], ['E82', 'K38c'], ['E35', 'L114'], ['E243', 'K248'], ['E243', 'N200c']]
loopC = [['E177', 'K148c'], ['D178', 'K148c'], ['E181', 'R133'], ['D185', 'I128']]
combs = regul + loopC

sims    = ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']
reps    = [1, 2, 3, 4]
chains  = ['A', 'B', 'C', 'D', 'E']
metrics = ['ecd_twist', 'beta_expansion', 'm2_m1_dist', 'nine_prime_dist', 'loopc', 'loopc2']

def stderr(array):
    return np.std(array) / np.sqrt(len(array))

matplotlib.rcParams.update({'font.size': 20})

for comb in combs:

    # PREPROCESSING

    fullResidueName = comb[0]           # E35
    fullTargetName  = comb[1]           # T158c
    residue = int(fullResidueName[1:])  # 35
    target  = fullTargetName[1:]        # 158c

    if target[-1] == 'p':
        chains2 = ['B', 'C', 'D', 'E', 'A']
        target = target[:-1]  # Remove trailing letter (158c -> 158).
    elif target[-1] == 'c':
        chains2 = ['E', 'A', 'B', 'C', 'D']
        target = target[:-1]  # Remove trailing letter (158c -> 158).
    else:
        chains2 = ['A', 'B', 'C', 'D', 'E']

    # GATHER DATA / full trajectories from 'bloom' directory.

    fullData   = makeSuperDict([sims, reps, chains, metrics, []])
    fullResult = makeSuperDict([sims, metrics, []])

    for sim in sims:
        for rep in reps:
            for chain in chains:
                fname = f'couple/all_all_{sim}_{rep}_{chain}_{chain}.txt'
                for idx in range(0, len(metrics)):
                    fullData[sim][rep][chain][metrics[idx]] = loadCol(fname, idx + 1, header=0)

    for sim in sims:
        for metric in metrics:
            for rep in reps:
                for chain in chains:
                    fullResult[sim][metric].append(np.mean(fullData[sim][rep][chain][metric]))

    # GATHER DATA / filtered data from the contact.

    filteredData   = makeSuperDict([sims, reps, chains, metrics, []])
    filteredResult = makeSuperDict([sims, metrics, []])

    for sim in sims:
        for rep in reps:
            for ii in range(0, len(chains)):
                fname = f'couple/{residue}_{target}_{sim}_{rep}_{chains[ii]}_{chains2[ii]}.txt'
                for idx in range(0, len(metrics)):
                        filteredData[sim][rep][chains[ii]][metrics[idx]] = loadCol(fname, idx + 1, header=0)

    for sim in sims:
        for metric in metrics:
            for rep in reps:
                for chain in chains:
                    # If a file is empty this will results in a np.mean being
                    # taken over an empty list, which will results in a np.NAN
                    # being returned. We need to remove these before we proceed.
                    array = filteredData[sim][rep][chain][metric]
                    if len(array) > 0:
                        filteredResult[sim][metric].append(np.mean(array))

    # MAKE THE BARPLOTS

    nameList = ['All', f'{fullResidueName}-{fullTargetName}']
    width    = 0.2
    x        = np.arange(len(nameList))

    for metric in metrics:

        meanList = {'4HFI_4': [], '4HFI_7': [], '6ZGD_4': [], '6ZGD_7': []}
        serrList = {'4HFI_4': [], '4HFI_7': [], '6ZGD_4': [], '6ZGD_7': []}
        for sim in sims:
            meanList[sim].append(np.mean(fullResult[sim][metric]))
            serrList[sim].append(stderr(fullResult[sim][metric]))
            meanList[sim].append(np.mean(filteredResult[sim][metric]))
            serrList[sim].append(stderr(filteredResult[sim][metric]))

        fig = plt.figure(figsize=(8, 6))
        ax  = fig.add_subplot()

        # 6ZGD_7
        mean4 = meanList['6ZGD_7']
        serr4 = serrList['6ZGD_7']
        ax.bar(     x - width * 1.5, mean4, width, color='#8856a7', label='closed, pH 7')
        ax.errorbar(x - width * 1.5, mean4, serr4, color='#8856a7', fmt='none', capsize=6, linewidth=2)

        # 6ZGD_4
        mean3 = meanList['6ZGD_4']
        serr3 = serrList['6ZGD_4']
        ax.bar(     x - width / 2.0, mean3, width, color='#9ebcda', label='closed, pH 4')
        ax.errorbar(x - width / 2.0, mean3, serr3, color='#9ebcda', fmt='none', capsize=6, linewidth=2)

        # 4HFI_7
        mean2 = meanList['4HFI_7']
        serr2 = serrList['4HFI_7']
        ax.bar(     x + width / 2.0, mean2, width, color='#8856a7', label='open, pH 7', edgecolor='w', lw=1, hatch='//')
        ax.errorbar(x + width / 2.0, mean2, serr2, color='#8856a7', fmt='none', capsize=6, linewidth=2)

        # 4HFI_4
        mean1 = meanList['4HFI_4']
        serr1 = serrList['4HFI_4']
        ax.bar(     x + width * 1.5, mean1, width, color='#9ebcda', label='open, pH 4', edgecolor='w', lw=1, hatch='\\\\')
        ax.errorbar(x + width * 1.5, mean1, serr1, color='#9ebcda', fmt='none', capsize=6, linewidth=2)

        ax.set_xticks(x, nameList)
        ax.legend(loc=1, prop={'size': 12})

        if metric == 'ecd_twist':
            plt.ylabel('ECD Twist (degrees)')
            plt.ylim(-12.5, -19)

        if metric == 'beta_expansion':
            plt.ylabel('Beta Expansion (Å)')
            plt.ylim(14, 17)

        if metric == 'm2_m1_dist':
            plt.ylabel('M2-M1(-) Distance (Å)')
            plt.ylim(11, 20)

        if metric == 'nine_prime_dist':
            plt.ylabel('9\' Radius (Å)')
            plt.ylim(2.5, 4)

        if metric == 'loopc':
            plt.ylabel('loopC-loopB Distance (Å)')
            plt.ylim(9, 16)

        if metric == 'loopc2':
            plt.ylabel('loopC-comp Distance (Å)')
            plt.ylim(13, 19)

        plt.legend()
        outname = f'couple/{residue}_{target}_{metric}.png'
        plt.tight_layout()
        plt.savefig(outname)
        os.system(f'convert {outname} -trim {outname}')
        plt.clf()
        plt.close()

    # MAKE THE HISTOGRAMS

    for metric in metrics:

        # Define metric-specific stuff

        if metric == 'ecd_twist':
            histRange = (-25, -5)
            xlabel = 'ECD Twist (degrees)'

        if metric == 'beta_expansion':
            histRange = (12, 18)
            xlabel = 'Beta Expansion (Å)'

        if metric == 'm2_m1_dist':
            histRange = (9, 21)
            xlabel = 'M2-M1(-) Distance (Å)'

        if metric == 'nine_prime_dist':
            histRange = (2, 5.5)
            xlabel = '9\' Radius (Å)'

        if metric == 'loopc':
            histRange = (8, 18)
            xlabel = 'loopC-loopB Distance (Å)'

        if metric == 'loopc2':
            histRange = (12, 20)
            xlabel = 'loopC-comp Distance (Å)'

        # Start building the figure

        fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(8, 8))

        for idx in [0, 1]:

            maxval = []
            subplt = axs[idx]

            for sim in sims:

                # Each figure contains four curves (one for each sim).

                histList = []

                for rep in reps:
                    for chain in chains:

                        # Get histogram values and corresponding bin_edges
                        if idx == 0:
                            x = fullData[sim][rep][chain][metric]
                        if idx == 1:
                            x = filteredData[sim][rep][chain][metric]

                        hist, bin_edges = np.histogram(x, density=True, bins=30, range=histRange)

                        # Another bug fix to prevent nan from np.mean([])
                        if 'nan' in [str(val) for val in hist]:
                            continue

                        histList.append(hist)

                # Compute mean and standard error.
                meanList = len(hist) * [0]
                errorList = len(hist) * [0]

                for kk in range(0, len(hist)):

                    # Create list of 20 values
                    temp = [0] * len(histList)  # 4*5=20
                    for ll in range(0, len(histList)):  # 4*5=20
                        temp[ll] = histList[ll][kk]

                    meanList[kk]  = np.mean(temp)
                    errorList[kk] = stderr(temp)

                # Plot mean and shaded region (standard error).
                A = []
                B = []
                for kk in range(0, len(meanList)):
                    A.append(meanList[kk] + errorList[kk])
                    B.append(meanList[kk] - errorList[kk])

                # The values in these dictionaries come from bloom/static.txt, which
                # in tern come from 4HFI_clean_crystal.pdb and 6ZGD_clean_cryoEM.pdb.
                # Margins in matplotlib are 1.05 * largest value = max(maxval).
                # This is to make sure that the vlines reach top of figure without rescaling.
                static4HFI = {'ecd_twist': -11.631, 'beta_expansion': 13.767, 'm2_m1_dist': 14.094, 'nine_prime_dist': 5.075, 'loopc': 0, 'loopc2': 0}
                static6ZGD = {'ecd_twist': -16.832, 'beta_expansion': 15.430, 'm2_m1_dist': 17.270, 'nine_prime_dist': 3.497, 'loopc': 0, 'loopc2': 0}

                if sim == '4HFI_4':
                    subplt.plot(bin_edges[1:], meanList, linewidth=2,   color='#9ebcda', label='open, pH 4', linestyle='--')
                    subplt.fill_between(bin_edges[1:], A, B, alpha=0.5, color='#9ebcda', hatch='\\\\', edgecolor='w', lw=1)

                if sim == '4HFI_7':
                    subplt.vlines(x=static4HFI[metric], ymin=0, ymax=10, color='black', label='open, crystal', linestyle='--')
                    subplt.plot(bin_edges[1:], meanList, linewidth=2,   color='#8856a7', label='open, pH 7', linestyle='--')
                    subplt.fill_between(bin_edges[1:], A, B, alpha=0.5, color='#8856a7', hatch='//',   edgecolor='w', lw=1)

                if sim == '6ZGD_4':
                    subplt.plot(bin_edges[1:], meanList, linewidth=2,   color='#9ebcda', label='closed, pH 4')
                    subplt.fill_between(bin_edges[1:], A, B, alpha=0.5, color='#9ebcda')

                if sim == '6ZGD_7':
                    subplt.vlines(x=static6ZGD[metric], ymin=0, ymax=10, color='black', label='closed, cryoEM')
                    subplt.plot(bin_edges[1:], meanList, linewidth=2,   color='#8856a7', label='closed, pH 7')
                    subplt.fill_between(bin_edges[1:], A, B, alpha=0.5, color='#8856a7')

                maxval.append(max(A))

            subplt.set_ylim(0, max(maxval) * 1.05)
            subplt.set_xlim(histRange)

            if idx == 0:
                subplt.set_ylabel('All')
                subplt.legend(loc='best', prop={'size': 14})

            if idx == 1:
                subplt.set_ylabel(f'{fullResidueName}-{fullTargetName}')
                subplt.set_xlabel(xlabel)

        outname = f'couple/{residue}_{target}_{metric}_hist.png'
        plt.tight_layout()
        plt.savefig(outname)
        os.system(f'convert {outname} -trim {outname}')
        plt.clf()
        plt.close()
