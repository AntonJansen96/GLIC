import math
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from science.parsing import loadCol
from science.utility import makeSuperDict

sims    = ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']
reps    = [1, 2, 3, 4]
chains  = ['A', 'B', 'C', 'D', 'E']
metrics = ['ecd_twist', 'ecd_spread', 'ecd_upper_spread', 'beta_expansion', 'm2_m1_dist', 'm2_radius', 'm1_kink', 'm1_kink_alt', 'nine_prime_dist', 'nine_prime_pore', 'minus_two_prime_dist', 'minus_two_prime_pore', 'c_loop']
fancy   = ['ECD twist', 'ECD spread', 'Upper ECD Spread', 'beta expansion', 'M1-M2 distance', 'M2 radius', 'M1 kink', 'M1 kink alter', '9\' distance', '9\' pore', '-2\' distance', '-2\' pore', 'loop C']

#! THIS PART CREATES THE .CSV FILE #############################################

df1 = pd.DataFrame()
df2 = pd.DataFrame()

superResult = makeSuperDict([sims, metrics, []])

for sim in sims:
    for rep in reps:
        for chain in chains:
            fname = 'bloom/{}_{}_{}.txt'.format(sim, rep, chain)
            for idx in range(0, len(metrics)):
                fname = f"/home/anton/GIT/GLIC/analysis/bloom/{sim}_{rep}_{chain}.txt"
                df1[f"{sim}_{rep}_{chain}_{metrics[idx]}"] = pd.Series(loadCol(fname, idx + 1, header=0))
                superResult[sim][metrics[idx]].append(np.mean(loadCol(fname, idx + 1, header=0)))

for sim in sims:
    for metric in metrics:
        df2[f"{sim}_{metric}"] = superResult[sim][metric]

df1.to_csv('domain_raw.csv')
df2.to_csv('domain_mean.csv')

#! THIS PART LOADS THE .CSV FILE ###############################################

df1 = pd.read_csv('domain_raw.csv')
df2 = pd.read_csv('domain_mean.csv')

#! THIS PART CREATES FIGURE 4 ##################################################

# Set global font size for figures.
matplotlib.rcParams.update({'font.size': 24})

def stderr(array):
    return np.std(array) / np.sqrt(len(array))

def doPlot1(actual, ylim, ylabel, fname):
    """This function makes the barplots (right side of figure in Bergh21 paper)."""

    # Initialize required data structures
    meanList = {}
    serrList = {}
    for key in sims:
        meanList[key] = []
        serrList[key] = []

    # Fill required data structures
    for sim in sims:
        for metric in actual:
            meanList[sim].append(np.mean(df2[f"{sim}_{metric}"]))
            serrList[sim].append(stderr(df2[f"{sim}_{metric}"]))

    # Get the fancy names and match for the elements of actual
    nameList = []
    for metric in actual:
        nameList.append(fancy[metrics.index(metric)])

    # Some plotting stuff
    width = 0.2
    x     = np.arange(len(nameList))
    fig   = plt.figure(figsize=(6, 6), dpi=300)
    ax    = fig.add_subplot()

    # 6ZGD_7
    mean4 = meanList['6ZGD_7']
    serr4 = serrList['6ZGD_7']
    ax.bar(     x - width * 1.5, mean4, width, color='#8856a7', label='closed, pH 7.0')
    ax.errorbar(x - width * 1.5, mean4, serr4, color='#8856a7', fmt='none', capsize=7, linewidth=3, markeredgewidth=2)

    # 6ZGD_4
    mean3 = meanList['6ZGD_4']
    serr3 = serrList['6ZGD_4']
    ax.bar(     x - width / 2.0, mean3, width, color='#9ebcda', label='closed, pH 4.0')
    ax.errorbar(x - width / 2.0, mean3, serr3, color='#9ebcda', fmt='none', capsize=7, linewidth=3, markeredgewidth=2)

    # 4HFI_7
    mean2 = meanList['4HFI_7']
    serr2 = serrList['4HFI_7']
    ax.bar(     x + width / 2.0, mean2, width, color='#8856a7', label='open, pH 7.0', edgecolor='w', lw=0, hatch='//', zorder=12)
    ax.errorbar(x + width / 2.0, mean2, serr2, color='#8856a7', fmt='none', capsize=7, linewidth=3, markeredgewidth=2, zorder=10)

    # 4HFI_4
    mean1 = meanList['4HFI_4']
    serr1 = serrList['4HFI_4']
    ax.bar(     x + width * 1.5, mean1, width, color='#9ebcda', label='open, pH 4.0', edgecolor='w', lw=0, hatch='\\\\', zorder=13)
    ax.errorbar(x + width * 1.5, mean1, serr1, color='#9ebcda', fmt='none', capsize=7, linewidth=3, markeredgewidth=2, zorder=11)

    ax.set_xticks([-0.3, -0.1, 0.1, 0.3], ['Closed\npH 7.0', 'Closed\npH 4.0', 'Open\npH 7.0', 'Open\npH 4.0'], rotation=0, size=16)

    plt.ylim(ylim)
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.savefig(fname)
    plt.clf()
    plt.close()

def doPlot2(actual, xlabel, fname):
    """This function makes the histograms (left side of figure in Bergh21 paper)."""

    plt.figure(figsize=(9, 6), dpi=300)
    maxval = []  # See below what this is used for.

    for sim in ['6ZGD_7', '6ZGD_4', '4HFI_7', '4HFI_4']:
        allData = []
        valuesList = []
        for rep in reps:
            for chain in chains:
                #! GET HISTOGRAM VALUES
                x = list(df1[f"{sim}_{rep}_{chain}_{metric}"])
                x = [val for val in x if not math.isnan(val)]  # Remove NaN values.
                allData += x
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

        # The values in these dictionaries come from bloom/static.txt, which
        # in tern come from 4HFI_clean_crystal.pdb and 6ZGD_clean_cryoEM.pdb.
        # Margins in matplotlib are 1.05 * largest value = max(maxval).
        # This is to make sure that the vlines reach top of figure without rescaling.
        static4HFI = {'ecd_twist': -11.631, 'ecd_spread': 23.796, 'ecd_upper_spread': 25.487, 'beta_expansion': 13.767, 'm2_m1_dist': 14.094, 'm2_radius': 21.583, 'm1_kink': 156.102, 'm1_kink_alt': 164.330, 'nine_prime_dist': 5.075, 'nine_prime_pore': 7.626, 'minus_two_prime_dist': 11.420, 'minus_two_prime_pore': 6.004, 'c_loop': 13.966}
        static6ZGD = {'ecd_twist': -16.832, 'ecd_spread': 23.548, 'ecd_upper_spread': 24.948, 'beta_expansion': 15.430, 'm2_m1_dist': 17.270, 'm2_radius': 19.425, 'm1_kink': 157.222, 'm1_kink_alt': 160.535, 'nine_prime_dist': 3.497, 'nine_prime_pore': 5.973, 'minus_two_prime_dist': 11.200, 'minus_two_prime_pore': 5.888, 'c_loop': 14.378}

        #* For some reason I don't fully understand, adding the following two
        #* lines fixed my bug.
        values, bins = np.histogram(allData, density=True, bins=20)
        plt.hist(values, bins, alpha=0.00)

        if sim == '4HFI_4':
            plt.plot(bins[1:], meanList, linewidth=2,   color='#9ebcda', label='open, pH 4.0', linestyle='--')
            plt.fill_between(bins[1:], A, B, alpha=0.5, color='#9ebcda', hatch='\\\\', edgecolor='w', lw=1)

        if sim == '4HFI_7':
            plt.vlines(x=static4HFI[actual], ymin=0, ymax=10, color='black', label='open, crystal', linestyle='--')
            plt.plot(bins[1:], meanList, linewidth=2,   color='#8856a7', label='open, pH 7.0', linestyle='--')
            plt.fill_between(bins[1:], A, B, alpha=0.5, color='#8856a7', hatch='//',   edgecolor='w', lw=1)

        if sim == '6ZGD_4':
            plt.plot(bins[1:], meanList, linewidth=2,   color='#9ebcda', label='closed, pH 4.0')
            plt.fill_between(bins[1:], A, B, alpha=0.5, color='#9ebcda')

        if sim == '6ZGD_7':
            plt.vlines(x=static6ZGD[actual], ymin=0, ymax=10, color='black', label='closed, cryoEM')
            plt.plot(bins[1:], meanList, linewidth=2,   color='#8856a7', label='closed, pH 7.0')
            plt.fill_between(bins[1:], A, B, alpha=0.5, color='#8856a7')

        maxval.append(max(A))

    plt.ylim(0, max(maxval) * 1.05)
    plt.xlabel(xlabel)
    plt.ylabel('Probability density')
    plt.tight_layout()
    plt.savefig(fname)
    plt.clf()
    plt.close()

doPlot1(['nine_prime_dist'],  [2.5, 3.75], '9\' Radius (Å)',     'fig1_right.png')
doPlot1(['m2_m1_dist'],       [14, 18], 'M2-M1(-) Distance (Å)', 'fig2_right.png')
doPlot1(['beta_expansion'],   [14, 16], 'Beta Expansion (Å)',    'fig3_right.png')
doPlot1(['ecd_twist'],       [-19, -13], 'ECD Twist (degrees)',  'fig4_right.png')

doPlot2('nine_prime_dist',  '9\' Radius (Å)',        'fig1_left.png')
doPlot2('m2_m1_dist',       'M2-M1(-) Distance (Å)', 'fig2_left.png')
doPlot2('beta_expansion',   'Beta Expansion (Å)',    'fig3_left.png')
doPlot2('ecd_twist',        'ECD Twist (degrees)',   'fig4_left.png')
