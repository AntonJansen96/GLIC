#!/usr/bin/env python3

import os
import MDAnalysis
import multiprocessing as mp
from science.utility import triplet2letter

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
pool = mp.Pool(processes=mp.cpu_count())
pool.starmap(task, items, chunksize=1)
