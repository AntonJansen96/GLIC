import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

sims    = ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']
reps    = [1, 2, 3, 4]
chains  = ['A', 'B', 'C', 'D', 'E']
targets = [13, 14, 26, 31, 32, 35, 49, 55, 67, 69, 75, 82, 86, 88, 91, 97, 104, 115, 122, 127, 136, 145, 147, 153, 154, 161, 163, 177, 178, 181, 185, 222, 235, 243, 272, 277, 282]
combs   = [('E26', ['V79p', 'N80p', 'V81p']), ('E35', ['L114', 'T158c', 'Na+']), ('E82', ['T36', 'K38c']), ('E177', ['K148c', 'R179']), ('E243', ['N200c', 'K248'])]

# THIS PART LOADS THE .CSV FILE ################################################

df1 = pd.read_csv('contacts.csv')
df2 = pd.read_csv('protonation.csv')

# THIS PART CREATES FIGURE 3 ###################################################

matplotlib.rcParams.update({'font.size': 26})

def stderr(array: list) -> float:
    return np.std(array) / np.sqrt(len(array))

for comb in combs:
    fullName = comb[0]
    include  = comb[1]
    target   = int(fullName[1:])

    # Initialize required data structures
    nameList = []
    meanList = {'4HFI_4': [], '4HFI_7': [], '6ZGD_4': [], '6ZGD_7': []}
    serrList = {'4HFI_4': [], '4HFI_7': [], '6ZGD_4': [], '6ZGD_7': []}

    # First set of bars is protonation:
    for sim in ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']:
        meanList[sim].append(np.mean(list(df2[f"{sim}_{target}"])))
        serrList[sim].append(stderr(list(df2[f"{sim}_{target}"])))

    # Gather occupancy data:
    for sim in ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']:
        for key in include:

            if key == 'Na+':
                continue

            meanList[sim].append(np.mean(list(df1[f"{target}_{sim}_{key}"])))
            serrList[sim].append(stderr(list(df1[f"{target}_{sim}_{key}"])))

            if sim == '4HFI_4':
                nameList.append(key)

    if 'Na+' in include:
        # We always want Na+ as the final set of bars
        nameList.append('Na+')
        for sim in ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']:
            try:
                meanList[sim].append(np.mean(list(df1[f"{target}_{sim}_Na+"])))
                serrList[sim].append(stderr(list(df1[f"{target}_{sim}_Na+"])))
            except KeyError:
                meanList[sim].append(0.0)
                serrList[sim].append(0.0)

    # We always want four sets of contact bars.

    for sim in sims:
        while len(meanList[sim]) != 4:
            meanList[sim].append(-1)
            serrList[sim].append(0.05)

    # This is a suggestion from Erik to make contact legend more clear.

    for idx in range(0, len(nameList)):
        nameList[idx] = nameList[idx]

    f, (a0, a1) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1, 3.7]}, figsize=(10, 5), dpi=300)

    width = 0.20

    mean4 = meanList['6ZGD_7']
    serr4 = serrList['6ZGD_7']
    mean3 = meanList['6ZGD_4']
    serr3 = serrList['6ZGD_4']
    mean2 = meanList['4HFI_7']
    serr2 = serrList['4HFI_7']
    mean1 = meanList['4HFI_4']
    serr1 = serrList['4HFI_4']

    # PROTONATION

    x = np.arange(1)

    a0.bar(     x - width * 1.5, mean4[0], width, color='#8856a7', label='closed, pH 7.0')
    a0.errorbar(x - width * 1.5, mean4[0], serr4[0], color='#8856a7', fmt='none', capsize=7, linewidth=3, markeredgewidth=2)

    a0.bar(     x - width / 2.0, mean3[0], width, color='#9ebcda', label='closed, pH 4.0')
    a0.errorbar(x - width / 2.0, mean3[0], serr3[0], color='#9ebcda', fmt='none', capsize=7, linewidth=3, markeredgewidth=2)

    a0.bar(     x + width / 2.0, mean2[0], width, color='#8856a7', label='open, pH 7.0', edgecolor='w', lw=0, hatch='//', zorder=12)
    a0.errorbar(x + width / 2.0, mean2[0], serr2[0], color='#8856a7', fmt='none', capsize=7, linewidth=3, markeredgewidth=2, zorder=11)

    a0.bar(     x + width * 1.5, mean1[0], width, color='#9ebcda', label='open, pH 4.0', edgecolor='w', lw=0, hatch='\\\\', zorder=13)
    a0.errorbar(x + width * 1.5, mean1[0], serr1[0], color='#9ebcda', fmt='none', capsize=7, linewidth=3, markeredgewidth=2, zorder=12)

    a0.set_xticks(x, [fullName])
    a0.set_ylim(0, 1)
    a0.set_ylabel('Protonation fraction')

    # CONTACTS

    x = np.arange(3)

    a1.bar(     x - width * 1.5, mean4[1:], width, color='#8856a7', label='closed, pH 7.0')
    a1.errorbar(x - width * 1.5, mean4[1:], serr4[1:], color='#8856a7', fmt='none', capsize=7, linewidth=3, markeredgewidth=2)

    a1.bar(     x - width / 2.0, mean3[1:], width, color='#9ebcda', label='closed, pH 4.0')
    a1.errorbar(x - width / 2.0, mean3[1:], serr3[1:], color='#9ebcda', fmt='none', capsize=7, linewidth=3, markeredgewidth=2)

    a1.bar(     x + width / 2.0, mean2[1:], width, color='#8856a7', label='open, pH 7.0', edgecolor='w', lw=0, hatch='//', zorder=12)
    a1.errorbar(x + width / 2.0, mean2[1:], serr2[1:], color='#8856a7', fmt='none', capsize=7, linewidth=3, markeredgewidth=2, zorder=11)

    a1.bar(     x + width * 1.5, mean1[1:], width, color='#9ebcda', label='open, pH 4.0', edgecolor='w', lw=0, hatch='\\\\', zorder=13)
    a1.errorbar(x + width * 1.5, mean1[1:], serr1[1:], color='#9ebcda', fmt='none', capsize=7, linewidth=3, markeredgewidth=2, zorder=12)

    if len(nameList) == 2:
        a1.set_xticks(x[0:-1], nameList)

    if len(nameList) == 3:
        a1.set_xticks(x, nameList)

    a1.set_ylim(0, 1)
    a1.set_ylabel('Occupancy')

    plt.tight_layout(pad=0.8)
    plt.savefig(f'{target}.png')
    plt.clf()
