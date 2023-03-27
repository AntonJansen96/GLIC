#!/usr/bin/env python3

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pickle
import os

from science.parsing import loadCol

# Set global font size for figures.
matplotlib.rcParams.update({'font.size': 26})

combs = [('E26', ['V79p', 'N80p', 'V81p']), ('E82', ['T36', 'K38c']), ('E177', ['K148c', 'R179']), ('D178', ['K148c', 'E177']), ('D32', ['Y119', 'R192']), ('E35', ['L114', 'T158c', 'Na+']), ('D122', ['R118', 'S123', 'D115']), ('E243', ['N200c', 'K248']), ('E222', ['Na+'])]

for comb in combs:
    print(comb)

    fullName = comb[0]
    include = comb[1]
    target = int(fullName[1:])

    # superData holds four dictionaries: 4HFI_4, 4HFI_7, 6ZGD_4, 6ZGD_7.
    # Each of these dictionaries contains the following key-value pairs:
    # key = the name of any residue encountered in the 20 sim_rep_target_chain files.
    # val = a list of occupancies containing reps x chains = 20 values.
    superData = {'4HFI_4': {}, '4HFI_7': {}, '6ZGD_4': {}, '6ZGD_7': {}}

    for sim in ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']:
        for rep in [1, 2, 3, 4]:
            for chain in ['A', 'B', 'C', 'D', 'E']:

                fileName = 'contacts/{}_{}_{}_{}.txt'.format(sim, rep, target, chain)
                list1 = loadCol(fileName, col=1, header=1)  # List of names ('T158')
                list2 = loadCol(fileName, col=2, header=1)  # List of chains ('A')
                list3 = loadCol(fileName, col=3, header=1)  # List of occs (0.951)

                # Change the chainIDs in list2 to either '', 'c', or 'p':
                for idx in range(0, len(list2)):
                    X = list2[idx]

                    # If the contact is in the same chain, remove the letter.
                    # The second conditional (chain = 'x') we add because 'x'
                    # is used as the 'default' chain for NA, CL, etc.
                    if (chain == X) or (X == 'x'):
                        list2[idx] = ''
                    # If we have this, the subunit containing the contact is complementary.
                    elif (chain == 'A' and X == 'E') or (chain == 'B' and X == 'A') or (chain == 'C' and X == 'B') or (chain == 'D' and X == 'C') or (chain == 'E' and X == 'D'):
                        list2[idx] = 'c'
                    # IF we have this, the subunit containing the contact is principal.
                    elif (chain == 'A' and X == 'B') or (chain == 'B' and X == 'C') or (chain == 'C' and X == 'D') or (chain == 'D' and X == 'E') or (chain == 'E' and X == 'A'):
                        list2[idx] = 'p'

                # Append the '', 'c', 'p' identifier to the names in list1:
                for idx in range(0, len(list1)):
                    list1[idx] += list2[idx]

                # Add the names and associated occupancies to superData:
                for idx in range(0, len(list1)):
                    # if there is already an entry for this name in superData
                    # simply add the occupancy to the existing list (value for said name).
                    if list1[idx] in superData[sim]:
                        superData[sim][list1[idx]].append(list3[idx])
                    # Else, create the key-value pair and initialize with [occ]
                    else:
                        superData[sim][list1[idx]] = [list3[idx]]

    # Fill any list not containing reps x chains = 20 elements up until 20
    # elements with additional zeros.
    for sim in ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']:
        for key in superData[sim]:
            if len(superData[sim][key]) != 20:
                superData[sim][key] += ((20 - len(superData[sim][key])) * [0.0])

    # Make sure all dictionaries contain all the same keys (initialize with lists
    # containing only zeros if necessary). This is to prevent KeyErrors later.
    sims = ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']
    for ii in range(0, len(sims)):
        for jj in range(0, len(sims)):
            if ii == jj:
                continue
            for key in superData[sims[ii]]:
                if key not in superData[sims[jj]]:
                    superData[sims[jj]][key] = 20 * [0.0]

    # Although the individual simulation dictionaries are sorted in descending
    # order, this does not mean that OVERALL (over the four sims) this should be
    # the order of the most occupied contact.

    totalList = {}                      # First create a list of key:value pairs
    for key in superData['4HFI_4']:     # where value is the sum of means of
                                        # occupancies for said key.
        total = 0
        for sim in ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']:
            total += np.mean(superData[sim][key])

        totalList[key] = float('{:.3f}'.format(total))

    # Sort dictionary key:value pairs by value (descending):
    totalList = dict(sorted(totalList.items(), key=lambda item: item[1], reverse=True))

    # Now use this totalList to also sort superData accordingly:
    # This is kind of hackish as I don't really understand lambdas in python
    # but it works.
    for sim in ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']:

        temp1 = sorted(superData[sim].items(), key=lambda kv: totalList[kv[0]], reverse=True)
        temp2 = {}

        for tuple in temp1:
            key        = tuple[0]
            value      = tuple[1]
            temp2[key] = value

        superData[sim] = temp2

    # print(superData['4HFI_4'])
    # The SuperData data structure should now be completed.

    #---------------------------------------------------------------------------

    protoMean = pickle.load(open('paperBarPlots.obj', 'rb'))
    protoData = {'4HFI_4': protoMean['4HFI_4'][target], '4HFI_7': protoMean['4HFI_7'][target], '6ZGD_4': protoMean['6ZGD_4'][target], '6ZGD_7': protoMean['6ZGD_7'][target]}

    #---------------------------------------------------------------------------

    # Here we do the actual plotting of the bar graphs.

    topNum = 5              # Base the top 'topNum'

    # Initialize required data structures
    nameList = []
    meanList = {'4HFI_4': [], '4HFI_7': [], '6ZGD_4': [], '6ZGD_7': []}
    serrList = {'4HFI_4': [], '4HFI_7': [], '6ZGD_4': [], '6ZGD_7': []}

    # First set of bars is protonation:
    for sim in ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']:
        meanList[sim].append(np.mean(protoData[sim]))
        serrList[sim].append(np.std(protoData[sim]) / np.sqrt(len(protoData[sim])))

    # Gather occupancy data:
    for sim in ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']:
        for key in include:

            if key == 'Na+':
                continue

            meanList[sim].append(np.mean(superData[sim][key]))
            serrList[sim].append(np.std(superData[sim][key]) / np.sqrt(len(superData[sim][key])))

            if sim == '4HFI_4':
                nameList.append(key)

    if 'Na+' in include:
        # We always want Na+ as the final set of bars
        nameList.append('Na+')
        for sim in ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']:
            try:
                meanList[sim].append(np.mean(superData[sim]['Na+']))
                serrList[sim].append(np.std(superData[sim]['Na+']) / len(superData[sim]['Na+']))
            except KeyError:
                meanList[sim].append(0.0)
                serrList[sim].append(0.0)

    #? #############################################################################

    #! We always want four sets of contact bars

    for sim in sims:
        while len(meanList[sim]) != 4:
            meanList[sim].append(-1)
            serrList[sim].append(0.05)

    #! This is a suggestion from Erik to make contact legend more clear.

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

    #! PROTONATION

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
    a0.set_ylabel('Protonation')

    #! CONTACTS

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

    # plt.legend(prop={'size': 16})
    plt.tight_layout(pad=0.8)
    plt.savefig(f'custom/{target}.png')
    os.system(f'convert custom/{target}.png -trim custom/{target}.png')
    plt.clf()
