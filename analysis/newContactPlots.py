#!/usr/bin/env python3

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pickle
from science.utility import makeSuperDict

matplotlib.rcParams.update({'font.size': 26})

# PARAMETERS #############################################################################

sims = ['6ZGD_7', '6ZGD_4', '4HFI_7', '4HFI_4']
residues = [13, 14, 26, 31, 32, 35, 49, 55, 67, 69, 75, 82, 86, 88, 91, 97, 104, 115, 122, 127, 136, 145, 147, 153, 154, 161, 163, 177, 178, 181, 185, 222, 235, 243, 272, 277, 282]

for residue in residues:

    #? GET TOP-4 MOST OCCUPIED CONTACTS ##################################################

    contacts = {}
    for sim in sims:
        with open(f"newcontacts/{sim}_{residue}.txt", 'r') as file:
            for line in file.read().splitlines():
                if line == '':
                    break

                identifier = line.split()[0]
                occupancy  = float(line.split()[1])

                if identifier not in contacts:
                    contacts[identifier] = occupancy
                else:
                    contacts[identifier] += occupancy

    # Get the identifiers of the top 4 most occupied residues
    top4 = [key for key in dict(sorted(contacts.items(), key=lambda item: item[1], reverse=True))][:4]
    # print(top4)  # debug.

    #? GET THE MEAN OCCUPANCIES FOR THE TOP-4 CONTACTS FOR ALL 4 SIMS ####################

    meanList = {'4HFI_4': [], '4HFI_7': [], '6ZGD_4': [], '6ZGD_7': []}
    serrList = {'4HFI_4': [], '4HFI_7': [], '6ZGD_4': [], '6ZGD_7': []}

    for sim in sims:
        with open(f"newcontacts/{sim}_{residue}.txt", 'r') as file:
            for identifier in top4:
                for line in file.read().splitlines():
                    if line.split()[0] == identifier:
                        meanList[sim].append(float(line.split()[1]))
                        file.seek(0)
                        break

    # print(meanList)  # debug.

    #? GET THE MEAN ATOM OCCUPANCIES FOR THE TOP-4 CONTACTS FOR ALL 4 SIMS ###############

    atomDict = makeSuperDict([sims, top4, {}])
    # print(atomDict)  # debug.
    for sim in sims:
        with open(f"newcontacts/{sim}_{residue}.txt", 'r') as file:
            for identifier in top4:
                for line in file.read().splitlines():
                    if len(line.split()) == 3 and line.split()[0] == identifier:
                        atom = line.split()[1]
                        occ = float(line.split()[2])
                        atomDict[sim][identifier][atom] = occ
                file.seek(0)

    # print(atomDict)  # debug.

    #? NORMALIZE TO 1 ####################################################################

    for sim in sims:
        for identifier in top4:
            total = 0
            for key in atomDict[sim][identifier]:
                total += atomDict[sim][identifier][key]
            factor = 1.0 / total
            for key in atomDict[sim][identifier]:
                atomDict[sim][identifier][key] *= factor

    # print(atomDict)  # debug

    #? SCALE TO GLOBAL OCCUPANCY #########################################################

    for sim in sims:
        for idx in range(0, len(top4)):
            globalMean = meanList[sim][idx]
            for key in atomDict[sim][top4[idx]]:
                atomDict[sim][top4[idx]][key] *= globalMean

    # print(atomDict)  # debug

    #? FORMAT TO WHAT THE INPUT LIST OF LISTS SHOULD LOOK LIKE ###########################

    #* 6ZGD_7
    maxAtoms = 0  # Maximum number of atoms any given residue contacts.
    for identifier in top4:
        L = len(atomDict['6ZGD_7'][identifier])
        if L > maxAtoms:
            maxAtoms = L

    mean4 = [[0, 0, 0, 0] for _ in range(maxAtoms)]
    mean4Labels = [['', '', '', ''] for _ in range(maxAtoms)]
    ii = 0
    for identifier in top4:
        jj = 0
        for key in atomDict['6ZGD_7'][identifier]:
            num = atomDict['6ZGD_7'][identifier][key]
            mean4[jj][ii] = num
            mean4Labels[jj][ii] = key
            jj += 1
        ii += 1
    # print(mean4)  # debug.

    #* 6ZGD_4
    maxAtoms = 0  # Maximum number of atoms any given residue contacts.
    for identifier in top4:
        L = len(atomDict['6ZGD_4'][identifier])
        if L > maxAtoms:
            maxAtoms = L

    mean3 = [[0, 0, 0, 0] for _ in range(maxAtoms)]
    mean3Labels = [['', '', '', ''] for _ in range(maxAtoms)]
    ii = 0
    for identifier in top4:
        jj = 0
        for key in atomDict['6ZGD_4'][identifier]:
            num = atomDict['6ZGD_4'][identifier][key]
            mean3[jj][ii] = num
            mean3Labels[jj][ii] = key
            jj += 1
        ii += 1
    # print(mean4)  # debug.

    #* 4HFI_7
    maxAtoms = 0  # Maximum number of atoms any given residue contacts.
    for identifier in top4:
        L = len(atomDict['4HFI_7'][identifier])
        if L > maxAtoms:
            maxAtoms = L

    mean2 = [[0, 0, 0, 0] for _ in range(maxAtoms)]
    mean2Labels = [['', '', '', ''] for _ in range(maxAtoms)]
    ii = 0
    for identifier in top4:
        jj = 0
        for key in atomDict['4HFI_7'][identifier]:
            num = atomDict['4HFI_7'][identifier][key]
            mean2[jj][ii] = num
            mean2Labels[jj][ii] = key
            jj += 1
        ii += 1
    # print(mean2)  # debug.

    #* 4HFI_4
    maxAtoms = 0  # Maximum number of atoms any given residue contacts.
    for identifier in top4:
        L = len(atomDict['4HFI_4'][identifier])
        if L > maxAtoms:
            maxAtoms = L

    mean1 = [[0, 0, 0, 0] for _ in range(maxAtoms)]
    mean1Labels = [['', '', '', ''] for _ in range(maxAtoms)]
    ii = 0
    for identifier in top4:
        jj = 0
        for key in atomDict['4HFI_4'][identifier]:
            num = atomDict['4HFI_4'][identifier][key]
            mean1[jj][ii] = num
            mean1Labels[jj][ii] = key
            jj += 1
        ii += 1
    # print(mean1)  # debug.

    #? TEST DATA #########################################################################

    # top4 = ['K248', 'R33', 'R177', 'R18']
    # atoms = ['C', 'N', 'OH1']

    # mean4 = [
    #     [0.13, 0.85, 0.25, 0.53],
    #     [0.66, 0.07, 0.25, 0.32],
    #     [0.00, 0.00, 0.13, 0.00]
    # ]

    # mean3 = [
    #     [0.13, 0.85, 0.25, 0.53],
    #     [0.66, 0.07, 0.25, 0.32],
    #     [0.00, 0.00, 0.13, 0.00]
    # ]

    # mean2 = [
    #     [0.13, 0.85, 0.25, 0.53],
    #     [0.66, 0.07, 0.25, 0.32],
    #     [0.00, 0.00, 0.13, 0.00]
    # ]

    # mean1 = [
    #     [0.13, 0.85, 0.25, 0.53],
    #     [0.66, 0.07, 0.25, 0.32],
    #     [0.00, 0.00, 0.13, 0.00]
    # ]

    f, (a0, a1) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1, 3.7]}, figsize=(18, 5), dpi=200)

    width = 0.22

    #! PLOTTING / PROTONATION ############################################################

    x = np.arange(1)

    protoMean = pickle.load(open('paperBarPlots.obj', 'rb'))
    protoData = protoData = {'4HFI_4': protoMean['4HFI_4'][residue], '4HFI_7': protoMean['4HFI_7'][residue], '6ZGD_4': protoMean['6ZGD_4'][residue], '6ZGD_7': protoMean['6ZGD_7'][residue]}
    protoMeanList = {'4HFI_4': [], '4HFI_7': [], '6ZGD_4': [], '6ZGD_7': []}
    protoSerrList = {'4HFI_4': [], '4HFI_7': [], '6ZGD_4': [], '6ZGD_7': []}

    for sim in sims:
        protoMeanList[sim].append(np.mean(protoData[sim]))
        protoSerrList[sim].append(np.std(protoData[sim]) / np.sqrt(len(protoData[sim])))

    m4 = protoMeanList['6ZGD_7']
    s4 = protoSerrList['6ZGD_7']
    m3 = protoMeanList['6ZGD_4']
    s3 = protoSerrList['6ZGD_4']
    m2 = protoMeanList['4HFI_7']
    s2 = protoSerrList['4HFI_7']
    m1 = protoMeanList['4HFI_4']
    s1 = protoSerrList['4HFI_4']

    a0.bar(     x - width * 1.5, m4, width, color='#8856a7', label='closed, pH 7.0')
    a0.errorbar(x - width * 1.5, m4, s4, color='#8856a7', fmt='none', capsize=7, linewidth=3, markeredgewidth=2)

    a0.bar(     x - width / 2.0, m3, width, color='#9ebcda', label='closed, pH 4.0')
    a0.errorbar(x - width / 2.0, m3, s3, color='#9ebcda', fmt='none', capsize=7, linewidth=3, markeredgewidth=2)

    a0.bar(     x + width / 2.0, m2, width, color='#8856a7', label='open, pH 7.0', edgecolor='w', lw=0, hatch='//', zorder=12)
    a0.errorbar(x + width / 2.0, m2, s2, color='#8856a7', fmt='none', capsize=7, linewidth=3, markeredgewidth=2, zorder=11)

    a0.bar(     x + width * 1.5, m1, width, color='#9ebcda', label='open, pH 4.0', edgecolor='w', lw=0, hatch='\\\\', zorder=13)
    a0.errorbar(x + width * 1.5, m1, s1, color='#9ebcda', fmt='none', capsize=7, linewidth=3, markeredgewidth=2, zorder=12)

    # u = MDAnalysis.Universe('../sims/4HFI_4/01/CA.pdb')
    # resname = u.select_atoms(f"chainID A and resid {residue}").residues.resnames[0]
    # resname = triplet2letter(resname) + str(residue)
    a0.set_xticks(x, [residue])

    a0.set_ylim(0, 1)
    a0.set_ylabel('Protonation')

    #! PLOTTING / CONTACTS ###############################################################

    x = np.arange(4)

    # 6ZGD_7
    ii = 0
    bottom = np.zeros(4)
    for weight in mean4:
        a1.bar(x - width * 1.5, weight, width, color='#8856a7', lw=1.5, edgecolor='white', bottom=bottom)
        for idx in range(0, len(x)):
            if weight[idx] > 0.05:  # This if-statement prevent plotting a label when the difference is very small.
                a1.text(x[idx] - width * 1.5 - 0.095, bottom[idx] + 0.01, mean4Labels[ii][idx], size=19, color='w')
        bottom += weight
        ii += 1

    # 6ZGD_4
    ii = 0
    bottom = np.zeros(4)
    for weight in mean3:
        a1.bar(x - width / 2.0, weight, width, color='#9ebcda', lw=1.5, edgecolor='white', bottom=bottom)
        for idx in range(0, len(x)):
            if weight[idx] > 0.05:  # This if-statement prevent plotting a label when the difference is very small.
                a1.text(x[idx] - width / 2.0 - 0.095, bottom[idx] + 0.01, mean3Labels[ii][idx], size=19, color='w')
        bottom += weight
        ii += 1

    # 4HFI_7
    ii = 0
    bottom = np.zeros(4)
    for weight in mean2:
        a1.bar(x + width / 2.0, weight, width, color='#8856a7', lw=1.5, edgecolor='white', hatch='//', bottom=bottom)
        for idx in range(0, len(x)):
            if weight[idx] > 0.05:  # This if-statement prevent plotting a label when the difference is very small.
                a1.text(x[idx] + width / 2.0 - 0.095, bottom[idx] + 0.01, mean2Labels[ii][idx], size=19, color='w')
        bottom += weight
        ii += 1

    # 4HFI_4
    ii = 0
    bottom = np.zeros(4)
    for weight in mean1:
        a1.bar(x + width * 1.5, weight, width, color='#9ebcda', lw=1.5, edgecolor='white', hatch='\\\\', bottom=bottom)
        for idx in range(0, len(x)):
            if weight[idx] > 0.05:  # This if-statement prevent plotting a label when the difference is very small.
                a1.text(x[idx] + width * 1.5 - 0.095, bottom[idx] + 0.01, mean1Labels[ii][idx], size=19, color='w')
        bottom += weight
        ii += 1

    a1.set_xticks(x, top4)
    a1.set_ylim(0, 1)
    a1.set_ylabel('Occupancy')

    plt.tight_layout(pad=0.4)
    plt.savefig(f'newcontacts/{residue}.png')
    plt.clf()
