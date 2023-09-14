import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from science.parsing import loadCol
from science.utility import makeSuperDict

sims    = ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']
reps    = [1, 2, 3, 4]
chains  = ['A', 'B', 'C', 'D', 'E']
targets = [13, 14, 26, 31, 32, 35, 49, 55, 67, 69, 75, 82, 86, 88, 91, 97, 104, 115, 122, 127, 136, 145, 147, 153, 154, 161, 163, 177, 178, 181, 185, 222, 235, 243, 272, 277, 282]
combs   = [('E26', ['V79p', 'N80p', 'V81p']), ('E35', ['L114', 'T158c', 'Na+']), ('E82', ['T36', 'K38c']), ('E177', ['K148c', 'R179']), ('E243', ['N200c', 'K248'])]

#! THIS PART CREATES THE .CSV FILE #############################################

df1 = pd.DataFrame()

for target in targets:

    # superData holds four dictionaries: 4HFI_4, 4HFI_7, 6ZGD_4, 6ZGD_7.
    # Each of these dictionaries contains the following key-value pairs:
    # key = the name of any residue encountered in the 20 sim_rep_target_chain files.
    # val = a list of occupancies containing reps x chains = 20 values.
    superData = makeSuperDict([sims, {}])

    for sim in sims:
        for rep in reps:
            for chain in chains:

                fileName = f"/home/anton/GIT/GLIC/analysis/contacts/{sim}_{rep}_{target}_{chain}.txt"
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
    for sim in sims:
        for key in superData[sim]:
            if len(superData[sim][key]) != 20:
                superData[sim][key] += ((20 - len(superData[sim][key])) * [0.0])

    # Make sure all dictionaries contain all the same keys (initialize with lists
    # containing only zeros if necessary). This is to prevent KeyErrors later.
    sims = sims
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
        for sim in sims:
            total += np.mean(superData[sim][key])

        totalList[key] = float('{:.3f}'.format(total))

    # Sort dictionary key:value pairs by value (descending):
    totalList = dict(sorted(totalList.items(), key=lambda item: item[1], reverse=True))

    # Now use this totalList to also sort superData accordingly:
    # This is kind of hackish as I don't really understand lambdas in python
    # but it works.
    for sim in sims:

        temp1 = sorted(superData[sim].items(), key=lambda kv: totalList[kv[0]], reverse=True)
        temp2 = {}

        for tuple in temp1:
            key        = tuple[0]
            value      = tuple[1]
            temp2[key] = value

        superData[sim] = temp2

    for sim in sims:
        for key in superData[sim]:
            df1[f"{target}_{sim}_{key}"] = superData[sim][key]

df1.to_csv('contacts.csv')

#! THIS PART LOADS THE .CSV FILE ###############################################

df1 = pd.read_csv('contacts.csv')
df2 = pd.read_csv('protonation.csv')

#! THIS PART CREATES FIGURE 3 ##################################################

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
    a0.set_ylabel('Protonation')

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
