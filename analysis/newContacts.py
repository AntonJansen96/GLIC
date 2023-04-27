#!/usr/bin/env python3

import os
import numpy as np
import pickle
import MDAnalysis
import multiprocessing as mp
import matplotlib
import matplotlib.pyplot as plt
from science.utility import triplet2letter

matplotlib.rcParams.update({'font.size': 26})

# PARAMETERS #############################################################################

residues = [13, 14, 26, 31, 32, 35, 49, 55, 67, 69, 75, 82, 86, 88, 91, 97, 104, 115, 122, 127, 136, 145, 147, 153, 154, 161, 163, 177, 178, 181, 185, 222, 235, 243, 272, 277, 282]
sims = ['6ZGD_7', '6ZGD_4', '4HFI_7', '4HFI_4']
reps = [1, 2, 3, 4]
chains = ['A', 'B', 'C', 'D', 'E']

carboxylAtoms = 'name OE1 OE2 OD1 OD2 NE2 ND1'
polarAtoms = 'name NZ NE2 OD2 OE1 OD1 NE OG1 OH ND2 OE2 OT1 NH2 O NH1 NE1 N ND1 OT2 OG'

cutoff = 4.0  # Angstroms.

##########################################################################################

# Task for multithreading.
def task(residue: int, sim: str) -> None:
    # Skip if analysis was already performed.
    if os.path.exists(f"newcontacts/{sim}_{residue}.txt"):
        return

    nFrames = 0
    atomDict = {}
    residueDict = {}

    for rep in reps:
        for chain in chains:

            path1 = f"../sims/{sim}/{rep:02d}/CA.pdb"
            path2 = f"../sims/{sim}/{rep:02d}/MD_conv.xtc"
            u = MDAnalysis.Universe(path1, path2)

            target = f"(resid {residue} and chainID {chain} and {carboxylAtoms})"
            consider = f"(chainID A B C D E and {polarAtoms} and not (resid {residue} and chainID {chain}))"

            for _ in u.trajectory:

                resSet = set()
                for atom in u.select_atoms(f"{consider} and around {cutoff} {target}"):

                    # Convert resid and chain information into a single string using
                    # complementary and principal identifiers.
                    X = atom.chainID
                    identifier = triplet2letter(atom.resname) + str(atom.resid)
                    # Case I: normal intrasubunit contact.
                    if chain == X:
                        pass
                    # Case II: subunit containing the contact is complementary.
                    elif (chain == 'A' and X == 'E') or (chain == 'B' and X == 'A') or (chain == 'C' and X == 'B') or (chain == 'D' and X == 'C') or (chain == 'E' and X == 'D'):
                        identifier += 'c'
                    # Case III: subunit containing the contact is principal.
                    elif (chain == 'A' and X == 'B') or (chain == 'B' and X == 'C') or (chain == 'C' and X == 'D') or (chain == 'D' and X == 'E') or (chain == 'E' and X == 'A'):
                        identifier += 'p'

                    atomTuple = (identifier, atom.name)

                    if atomTuple not in atomDict:
                        atomDict[atomTuple] = 1
                    else:
                        atomDict[atomTuple] += 1

                    resSet.add(identifier)

                for element in resSet:
                    if element not in residueDict:
                        residueDict[element] = 1
                    else:
                        residueDict[element] += 1

                nFrames += 1

    # Make ATOM occupancies fractional by dividing by #hits residue received in total.
    for key1 in residueDict:
        for key2 in atomDict:
            if key1 == key2[0]:
                atomDict[key2] /= float(residueDict[key1])

    # Make RESIDUE occupancies fractional by dividing #hits by total #frames.
    for key in residueDict:
        residueDict[key] /= float(nFrames)

    # Sort the the residueDict by value
    residueDict = dict(sorted(residueDict.items(), key=lambda item: item[1], reverse=True))

    # Write residue occupancies to file.
    with open(f"newcontacts/{sim}_{residue}.txt", 'w') as file:
        for key in residueDict:
            file.write(f"{key:<5s} {residueDict[key]:.3f}\n")
        # Write atom occupancies to file.
        file.write('\n')
        for key1 in residueDict:
            for key2 in atomDict:
                if key1 == key2[0]:
                    file.write(f"{key1:<5s} {key2[1]:<4s} {atomDict[key2]:.3f}\n")

# Gather items:
items = []
for residue in residues:
    for sim in sims:
        items.append((residue, sim))

# Run multithreaded
# pool = mp.Pool(processes=mp.cpu_count())
# pool.starmap(task, items, chunksize=1)

##########################################################################################

protoMean = pickle.load(open('paperBarPlots.obj', 'rb'))

for residue in [35]:

    protoData = {'4HFI_4': protoMean['4HFI_4'][residue], '4HFI_7': protoMean['4HFI_7'][residue], '6ZGD_4': protoMean['6ZGD_4'][residue], '6ZGD_7': protoMean['6ZGD_7'][residue]}

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

    print(top4)

    ######################################################################################
    # FILL meanList. First element is protonation, last 4 elements are contacts.

    meanList = {'4HFI_4': [], '4HFI_7': [], '6ZGD_4': [], '6ZGD_7': []}
    serrList = {'4HFI_4': [], '4HFI_7': [], '6ZGD_4': [], '6ZGD_7': []}

    for sim in sims:
        meanList[sim].append(np.mean(protoData[sim]))
        serrList[sim].append(np.std(protoData[sim]) / np.sqrt(len(protoData[sim])))

    for sim in sims:
        with open(f"newcontacts/{sim}_{residue}.txt", 'r') as file:
            for identifier in top4:
                for line in file.read().splitlines():
                    if line.split()[0] == identifier:
                        meanList[sim].append(float(line.split()[1]))
                        file.seek(0)
                        break

    ######################################################################################
    # PLOTTING

    f, (a0, a1) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1, 3.7]}, figsize=(11, 5), dpi=200)

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

    # u = MDAnalysis.Universe('../sims/4HFI_4/01/CA.pdb')
    # resname = u.select_atoms(f"chainID A and resid {residue}").residues.resnames[0]
    # resname = triplet2letter(resname) + str(residue)
    a0.set_xticks(x, [residue])

    a0.set_ylim(0, 1)
    a0.set_ylabel('Protonation')

    #! CONTACTS

    x = np.arange(4)

    a1.bar(     x - width * 1.5, mean4[1:], width, color='#8856a7', label='closed, pH 7.0')
    # a1.errorbar(x - width * 1.5, mean4[1:], serr4[1:], color='#8856a7', fmt='none', capsize=7, linewidth=3, markeredgewidth=2)
    a1.text((0 - width * 1.5 - 0.065), 0.02, "test1", rotation=90, size=18)
    a1.text((1 - width * 1.5 - 0.065), 0.02, "test2", rotation=90, size=18)
    a1.text((2 - width * 1.5 - 0.065), 0.02, "test3", rotation=90, size=18, zorder=20)
    a1.text((3 - width * 1.5 - 0.065), 0.02, "test4", rotation=90, size=18, zorder=20)

    a1.bar(     x - width / 2.0, mean3[1:], width, color='#9ebcda', label='closed, pH 4.0')
    # a1.errorbar(x - width / 2.0, mean3[1:], serr3[1:], color='#9ebcda', fmt='none', capsize=7, linewidth=3, markeredgewidth=2)
    a1.text((0 - width / 2.0 - 0.065), 0.02, "test5", rotation=90, size=18)
    a1.text((1 - width / 2.0 - 0.065), 0.02, "test6", rotation=90, size=18)
    a1.text((2 - width / 2.0 - 0.065), 0.02, "test7", rotation=90, size=18, zorder=20)
    a1.text((3 - width / 2.0 - 0.065), 0.02, "test8", rotation=90, size=18, zorder=20)

    a1.bar(     x + width / 2.0, mean2[1:], width, color='#8856a7', label='open, pH 7.0', edgecolor='w', lw=0, hatch='//', zorder=12)
    # a1.errorbar(x + width / 2.0, mean2[1:], serr2[1:], color='#8856a7', fmt='none', capsize=7, linewidth=3, markeredgewidth=2, zorder=11)
    a1.text((0 + width / 2.0 - 0.065), 0.02, "test9", rotation=90, size=18)
    a1.text((1 + width / 2.0 - 0.065), 0.02, "test10", rotation=90, size=18)
    a1.text((2 + width / 2.0 - 0.065), 0.02, "test11", rotation=90, size=18, zorder=20)
    a1.text((3 + width / 2.0 - 0.065), 0.02, "test12", rotation=90, size=18, zorder=20)

    a1.bar(     x + width * 1.5, mean1[1:], width, color='#9ebcda', label='open, pH 4.0', edgecolor='w', lw=0, hatch='\\\\', zorder=13)
    # a1.errorbar(x + width * 1.5, mean1[1:], serr1[1:], color='#9ebcda', fmt='none', capsize=7, linewidth=3, markeredgewidth=2, zorder=12)
    a1.text((0 + width * 1.5 - 0.065), 0.02, "test13", rotation=90, size=18)
    a1.text((1 + width * 1.5 - 0.065), 0.02, "test14", rotation=90, size=18)
    a1.text((2 + width * 1.5 - 0.065), 0.02, "test15", rotation=90, size=18, zorder=20)
    a1.text((3 + width * 1.5 - 0.065), 0.02, "test16", rotation=90, size=18, zorder=20)

    a1.set_xticks(x, top4)
    a1.set_ylim(0, 1)
    a1.set_ylabel('Occupancy')

    plt.tight_layout(pad=0.8)
    plt.savefig(f'newcontacts/{residue}.png')
    plt.clf()
