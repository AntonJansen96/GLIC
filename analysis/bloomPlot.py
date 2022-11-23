#!/bin/python3

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from science.utility import makeSuperDict
from science.parsing import loadCol

# Set global font size for figures.
matplotlib.rcParams.update({'font.size': 14})


def stderr(array):
    return np.std(array) / np.sqrt(len(array))


sims    = ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']
reps    = [1, 2, 3, 4]
chains  = ['A', 'B', 'C', 'D', 'E']
metrics = ['ecd_twist', 'ecd_spread', 'ecd_upper_spread', 'beta_expansion', 'm2_m1_dist', 'm2_radius', 'm1_kink', 'm1_kink_alt', 'nine_prime_dist', 'nine_prime_pore', 'minus_two_prime_dist', 'minus_two_prime_pore', 'c_loop']
fancy   = ['ECD twist', 'ECD spread', 'Upper ECD Spread', 'beta expansion', 'M1-M2 distance', 'M2 radius', 'M1 kink', 'M1 kink alter', '9\' distance', '9\' pore', '-2\' distance', '-2\' pore', 'loop C']

superData = makeSuperDict([sims, reps, chains, metrics, []])
superResult = makeSuperDict([sims, metrics, []])

# Load the data.
for sim in sims:
    for rep in reps:
        for chain in chains:
            fmetric = 'bloom/{}_{}_{}.txt'.format(sim, rep, chain)
            for idx in range(0, len(metrics)):
                superData[sim][rep][chain][metrics[idx]] = loadCol(fmetric, idx + 1, header=0)

# Get the results.
for sim in sims:
    for metric in metrics:
        for rep in reps:
            for chain in chains:
                superResult[sim][metric].append(np.mean(superData[sim][rep][chain][metric]))

#? PLOTTING PART OF SCRIPT


def doPlot1(actual, ylim, ylabel, fname):
    """This function makes the barplots (right side of figure in Bergh21 paper).
    """

    # Initialize required data structures
    meanList = {}
    serrList = {}
    for key in sims:
        meanList[key] = []
        serrList[key] = []

    # Fill required data structures
    for sim in sims:
        for metric in actual:
            meanList[sim].append(np.mean(superResult[sim][metric]))
            serrList[sim].append(stderr(superResult[sim][metric]))

    # Get the fancy names and match for the elements of actual
    nameList = []
    for metric in actual:
        nameList.append(fancy[metrics.index(metric)])

    # Some plotting stuff
    width = 0.2
    x     = np.arange(len(nameList))
    fig   = plt.figure(figsize=(6, 6))
    ax    = fig.add_subplot()

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

    # ax.set_xticks(x, nameList, rotation=0, color='white')
    ax.set_xticks([-0.3, -0.1, 0.1, 0.3], ['Closed\npH 7', 'Closed\npH 4', 'Open\npH 7', 'Open\npH 4'], rotation=0)
    # ax.legend(loc='best', prop={'size': 12})

    # plt.title(title)
    plt.ylim(ylim)
    plt.ylabel(ylabel, size=19)
    plt.tight_layout()
    plt.savefig('bloom/{}'.format(fname))
    plt.clf()
    plt.close()


def doPlot2(actual, xlabel, fname):
    """This function makes the histograms (left side of figure in Bergh21 paper).
    """

    plt.figure(figsize=(9, 6))

    for sim in ['6ZGD_7', '6ZGD_4', '4HFI_7', '4HFI_4']:
        valuesList = []
        for rep in reps:
            for chain in chains:

                #! GET HISTOGRAM VALUES
                x = superData[sim][rep][chain][actual]
                values, bins = np.histogram(x, density=True, bins=20)
                valuesList.append(values)

        #! COMPUTE MEAN AND STANDARD ERROR
        meanList  = len(values) * [0]  # 200, to hold mean for each bin
        errorList = len(values) * [0]  # 200, to hold erro for each bin

        for kk in range(0, len(values)):  # 200

            # Create list of 20 values
            temp = [0] * len(valuesList)  # 4*5=20
            for ll in range(0, len(valuesList)):  # 4*5=20
                temp[ll] = valuesList[ll][kk]

            meanList[kk]  = np.mean(temp)
            errorList[kk] = stderr(temp)

        #! PLOT MEAN AND SHADED REGION (STANDARD ERROR)
        A = []
        B = []
        for kk in range(0, len(meanList)):
            A.append(meanList[kk] + errorList[kk])
            B.append(meanList[kk] - errorList[kk])

        if sim == '4HFI_4':
            plt.plot(bins[1:], meanList, linewidth=2,   color='#9ebcda', label='open, pH 4', linestyle='--')
            plt.fill_between(bins[1:], A, B, alpha=0.5, color='#9ebcda', hatch='\\\\', edgecolor='w', lw=1)

        if sim == '4HFI_7':
            plt.plot(bins[1:], meanList, linewidth=2,   color='#8856a7', label='open, pH 7', linestyle='--')
            plt.fill_between(bins[1:], A, B, alpha=0.5, color='#8856a7', hatch='//',   edgecolor='w', lw=1)

        if sim == '6ZGD_4':
            plt.plot(bins[1:], meanList, linewidth=2,   color='#9ebcda', label='closed, pH 4')
            plt.fill_between(bins[1:], A, B, alpha=0.5, color='#9ebcda')

        if sim == '6ZGD_7':
            plt.plot(bins[1:], meanList, linewidth=2,   color='#8856a7', label='closed, pH 7')
            plt.fill_between(bins[1:], A, B, alpha=0.5, color='#8856a7')

    plt.xlabel(xlabel, size=19)
    plt.ylabel('Occupancy', size=19)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig('bloom/{}'.format(fname))
    plt.clf()
    plt.close()


doPlot1(['m2_radius'],        [19, 21], 'M2 Spread (Å)',         'fig1_right.png')  # good
doPlot1(['nine_prime_pore'],  [6,   7], '9\' Radius (Å)',        'fig2_right.png')  # good
doPlot1(['m2_m1_dist'],       [14, 18], 'M2-M1(-) Distance (Å)', 'fig3_right.png')  # good
doPlot1(['beta_expansion'],   [14, 16], 'Beta Expansion (Å)',    'fig4_right.png')  # good
doPlot1(['ecd_spread'],       [23, 26], 'ECD Spread (Å)',        'fig5_right.png')  # good
doPlot1(['ecd_upper_spread'], [25, 28], 'Upper ECD Spread (Å)',  'fig6_right.png')  # good

doPlot2('m2_radius',        'M2 Spread (Å)',         'fig1_left.png')  # good
doPlot2('nine_prime_pore',  '9\' Radius (Å)',        'fig2_left.png')  # good
doPlot2('m2_m1_dist',       'M2-M1(-) Distance (Å)', 'fig3_left.png')  # good
doPlot2('beta_expansion',   'Beta Expansion (Å)',    'fig4_left.png')  # good
doPlot2('ecd_spread',       'ECD Spread (Å)',        'fig5_left.png')  # good
doPlot2('ecd_upper_spread', 'Upper ECD Spread (Å)',  'fig6_left.png')  # good

#? Sanity check: this should be the same as the plot from the original script.
# for chain in chains:
#     x = superData['4HFI_4'][1][chain]['m2_radius']
#     t = range(0, len(x))
#     plt.plot(t, x, label=chain)
# plt.legend()
# plt.xlabel('Time (ns)')
# plt.ylabel('Distance (A??)')
# plt.title('4HFI_4 1 chain A m2_radius')
# plt.tight_layout()
# plt.savefig('bloom/m2_radius_comp.png')
