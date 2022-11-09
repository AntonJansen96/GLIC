#!/bin/python3

import matplotlib
import science as sc

# Set global font size for figures.
matplotlib.rcParams.update({'font.size': 14})


# parameters
topContacts = 5

# superData holds four dictionaries: 4HFI_4, 4HFI_7, 6ZGD_4, 6ZGD_7.
# Each of these four dictionaries contains the follwing key-value pairs:
# key = the name of any residue encountered in the 20 list.txt files.
# val = a list of occupancies.
superData = {'4HFI_4': {}, '4HFI_7': {}, '6ZGD_4': {}, '6ZGD_7': {}}

for target in [35]:
    for sim in ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']:
        for rep in [1, 2, 3, 4]:
            for chain in ['A', 'B', 'C', 'D', 'E']:

                fileName = 'contacts/{}_{}_{}_{}.txt'.format(sim, rep, target, chain)
                list1 = sc.parsing.loadCol(fileName, col=1, header=1)  # List of names ('T158')
                list2 = sc.parsing.loadCol(fileName, col=2, header=1)  # List of chains ('A')
                list3 = sc.parsing.loadCol(fileName, col=3, header=1)  # List of occs (0.951)

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
                    else:
                        print("how the fuck did we end up here?")

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
