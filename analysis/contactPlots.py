#!/bin/python3

import matplotlib
import numpy as np

from science.cphmd import getLambdaFileIndices
from science.cphmd import protonation
from science.parsing import loadCol
from science.parsing import loadxvg


# Set global font size for figures.
matplotlib.rcParams.update({'font.size': 14})

for target in [35]:

    # superData holds four dictionaries: 4HFI_4, 4HFI_7, 6ZGD_4, 6ZGD_7.
    # Each of these four dictionaries contains the follwing key-value pairs:
    # key = the name of any residue encountered in the 20 list.txt files.
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
                    if chain == X:
                        list2[idx] = ''
                    # If we have this, the subunit containing the contact is complementary.
                    elif (chain == 'A' and X == 'E') or (chain == 'B' and X == 'A') or (chain == 'C' and X == 'B') or (chain == 'D' and X == 'C') or (chain == 'E' and X == 'D'):
                        list2[idx] = 'c'
                    # IF we have this, the subunit containing the contact is principal.
                    elif (chain == 'A' and X == 'B') or (chain == 'B' and X == 'C') or (chain == 'C' and X == 'D') or (chain == 'D' and X == 'E') or (chain == 'E' and X == 'A'):
                        list2[idx] = 'p'

                # Add the '', 'c', 'p' identifier to the names in list1:
                for idx in range(0, len(list1)):
                    list1[idx] += list2[idx]

                # Add the names and associated occupancies to superData.
                for idx in range(0, len(list1)):
                    # if there is already an entry for this name in superData
                    # simply add the occupancy to the existing list (value for said name).
                    if list1[idx] in superData[sim]:
                        superData[sim][list1[idx]].append(list3[idx])
                    # Else, create the key-value pair:
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

    # print(superData['4HFI_4'])
    # The SuperData data structure should now be completed.

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
