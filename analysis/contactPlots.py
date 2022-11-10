#!/bin/python3

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import MDAnalysis

from science.cphmd import getLambdaFileIndices
from science.cphmd import protonation
from science.parsing import loadCol
from science.parsing import loadxvg
from science.utility import triplet2letter

# Set global font size for figures.
matplotlib.rcParams.update({'font.size': 14})

# 145, 161, 185, 282 are bad (have fewer than 5 inter-subunit contacts)

for target in [13, 14, 26, 31, 32, 35, 49, 55, 67, 69, 75, 82, 86, 88, 91, 97, 104, 115, 122, 127, 136, 145, 147, 153, 154, 161, 163, 177, 178, 181, 185, 222, 235, 243, 272, 277, 282]:

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

    # protoData is a dictionary that holds four lists. Each list contains 20
    # protonation values, where each value is the average protonation of a rep x chain.
    protoData = {'4HFI_4': [], '4HFI_7': [], '6ZGD_4': [], '6ZGD_7': []}
    lambdaIndices = getLambdaFileIndices('../sims/4HFI_4/01/CA.pdb', target)

    for sim in ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']:
        for rep in [1, 2, 3, 4]:
            for idx in lambdaIndices:

                x = loadxvg('../sims/{}/{:02d}/cphmd-coord-{}.xvg'.format(sim, rep, idx), dt=1000, b=0)[1]
                protoData[sim].append(protonation(x))

    # print(protoData)
    # The protoData data structure should now be completed.

    #---------------------------------------------------------------------------

    # Print our data structures to files for checking/debugging purposes.
    with open('contacts/{}.dat'.format(target), 'w') as file:
        file.write('#      | 4HFI_4                    | 4HFI_7                     | 6ZGD_4                     | 6ZGD_4\n')
        file.write('# name | |pro| serr   |occ| serr   | |pro| serr   |occ|  serr   | |pro| serr   |occ|  serr   | |pro| serr   |occ|  serr\n')
        for name in superData['4HFI_4']:
            file.write('{:<6s} | {:.3f} {:.4f} {:.3f} {:.4f} | {:.3f} {:.4f} {:.4f} {:.4f} | {:.3f} {:.4f} {:.4f} {:.4f} | {:.3f} {:.4f} {:.4f} {:.4f}\n'.format(
                name,
                np.mean(protoData['4HFI_4']),
                np.std( protoData['4HFI_4']) / np.sqrt(len(protoData['4HFI_4'])),
                np.mean(superData['4HFI_4'][name]),
                np.std( superData['4HFI_4'][name]) / np.sqrt(len(superData['4HFI_4'][name])),
                np.mean(protoData['4HFI_7']),
                np.std( protoData['4HFI_7']) / np.sqrt(len(protoData['4HFI_7'])),
                np.mean(superData['4HFI_7'][name]),
                np.std( superData['4HFI_7'][name]) / np.sqrt(len(superData['4HFI_7'][name])),
                np.mean(protoData['6ZGD_4']),
                np.std( protoData['6ZGD_4']) / np.sqrt(len(protoData['6ZGD_4'])),
                np.mean(superData['6ZGD_4'][name]),
                np.std( superData['6ZGD_4'][name]) / np.sqrt(len(superData['6ZGD_4'][name])),
                np.mean(protoData['6ZGD_7']),
                np.std( protoData['6ZGD_7']) / np.sqrt(len(protoData['6ZGD_7'])),
                np.mean(superData['6ZGD_7'][name]),
                np.std( superData['6ZGD_7'][name]) / np.sqrt(len(superData['6ZGD_7'][name]))
            ))

    #---------------------------------------------------------------------------

    # Here we do the actual plotting of the bar graphs.

    topNum = 5              # Base the top 'topNum'
    inWhich = '4HFI_4'      # contacts on sim 'inWhich'...

    # Initialize required data structures
    nameList = []
    meanList = {'4HFI_4': [], '4HFI_7': [], '6ZGD_4': [], '6ZGD_7': []}
    serrList = {'4HFI_4': [], '4HFI_7': [], '6ZGD_4': [], '6ZGD_7': []}

    # First set of bars is protonation:
    nameList.append('protonation')
    for sim in ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']:
        meanList[sim].append(np.mean(protoData[sim]))
        serrList[sim].append(np.std(protoData[sim]) / np.sqrt(len(protoData[sim])))

    # Gather occupancy data:
    count = 0
    for sim in ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']:
        for key in superData[inWhich]:

            if key == 'Na+':  # This is to prevent Na+ from showing up in the
                continue      # middle of the bar plot if Na+ is in the top 5.

            if sim == inWhich:
                nameList.append(key)

            meanList[sim].append(np.mean(superData[sim][key]))
            serrList[sim].append(np.std(superData[sim][key]) / np.sqrt(len(superData[sim][key])))

            count += 1
            if count == topNum:
                count = 0
                break

    # We always want Na+ as the final set of bars
    nameList.append('Na+')
    for sim in ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']:
        try:
            meanList[sim].append(np.mean(superData[sim]['Na+']))
            serrList[sim].append(np.std(superData[sim]['Na+']) / len(superData[sim]['Na+']))
        except KeyError:
            meanList[sim].append(0.0)
            serrList[sim].append(0.0)

    # Some plotting stuff
    width = 0.2
    x     = np.arange(len(nameList))
    fig   = plt.figure(figsize=(10, 6))
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

    ax.set_xticks(x, nameList)
    ax.legend(loc=1, prop={'size': 12})

    # Vertical line to separate protonation bars from rest
    plt.axvline(x=0.5, ymin=0, ymax=1.1, color='black', lw=2, linestyle=':')

    # Do the title
    u = MDAnalysis.Universe('../sims/4HFI_4/01/CA.pdb')
    triplet = u.select_atoms('chainID A and resid {}'.format(target)).residues.resnames[0]
    letter = triplet2letter(triplet)
    plt.title('Residue {}{} (top {} contacts overall)'.format(letter, target, topNum))

    plt.ylim(0, 1.1)
    plt.ylabel('Protonation / contact occupancy')
    plt.tight_layout()
    plt.savefig('contacts/{}.png'.format(target))
    plt.clf()
    plt.close()

    #---------------------------------------------------------------------------

    # Make another bar plot, but now for just the INTER-subunit contacts
    # (so either 'c' or 'p' relative to target).

    # Initialize required data structures
    nameList = []
    meanList = {'4HFI_4': [], '4HFI_7': [], '6ZGD_4': [], '6ZGD_7': []}
    serrList = {'4HFI_4': [], '4HFI_7': [], '6ZGD_4': [], '6ZGD_7': []}

    # First set of bars is protonation:
    nameList.append('protonation')
    for sim in ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']:
        meanList[sim].append(np.mean(protoData[sim]))
        serrList[sim].append(np.std(protoData[sim]) / np.sqrt(len(protoData[sim])))

    # Gather occupancy data:
    count = 0
    for sim in ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']:
        for key in superData[inWhich]:

            # This if-statement makes sure we only get inter-subunit contacts.
            if key[-1] not in ['c', 'p']:
                continue

            if key == 'Na+':  # This is to prevent Na+ from showing up in the
                continue      # middle of the bar plot if Na+ is in the top 5.

            if sim == inWhich:
                nameList.append(key)

            meanList[sim].append(np.mean(superData[sim][key]))
            serrList[sim].append(np.std(superData[sim][key]) / np.sqrt(len(superData[sim][key])))

            count += 1
            if count == topNum:
                count = 0
                break

    # We always want Na+ as the final set of bars
    nameList.append('Na+')
    for sim in ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']:
        try:
            meanList[sim].append(np.mean(superData[sim]['Na+']))
            serrList[sim].append(np.std(superData[sim]['Na+']) / len(superData[sim]['Na+']))
        except KeyError:
            meanList[sim].append(0.0)
            serrList[sim].append(0.0)

    # Some plotting stuff
    width = 0.2
    x     = np.arange(len(nameList))
    fig   = plt.figure(figsize=(10, 6))
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

    ax.set_xticks(x, nameList)
    ax.legend(loc=1, prop={'size': 12})

    # Vertical line to separate protonation bars from rest
    plt.axvline(x=0.5, ymin=0, ymax=1.1, color='black', lw=2, linestyle=':')

    # Do the title
    plt.title('Residue {}{} (top {} inter-subunit contacts)'.format(letter, target, topNum))

    plt.ylim(0, 1.1)
    plt.ylabel('Protonation / contact occupancy')
    plt.tight_layout()
    plt.savefig('contacts/{}_inter.png'.format(target))
    plt.clf()
    plt.close()
