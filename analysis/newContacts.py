#!/usr/bin/env python3

import MDAnalysis
import os
from tqdm import tqdm

# PARAMETERS #############################################################################

sims = ['6ZGD_7']
reps = [1, 2, 3, 4]
residues = [35]
chains = ['A']

carboxylAtoms = 'name OE1 OE2 OD1 OD2 NE2 ND1'
polarAtoms = 'name NZ NE2 OD2 OE1 OD1 NE OG1 OH ND2 OE2 OT1 NH2 O NH1 NE1 N ND1 OT2 OG'

cutoff = 4.0  # Angstroms.
ignoreNeighbours = 0  # implement this feature later.

##########################################################################################

if not os.path.exists('newcontacts'):
    os.system('mkdir newcontacts')

for residue in residues:
    for sim in sims:
        for chain in chains:

            nFrames = 0
            atomDict = {}
            residueDict = {}

            for rep in reps:

                path1 = f"../sims/{sim}/{rep:02d}/CA.pdb"
                path2 = f"../sims/{sim}/{rep:02d}/MD_conv.xtc"
                u = MDAnalysis.Universe(path1, path2)

                sel1str = f"(resid {residue} and chainID {chain} and {carboxylAtoms})"
                sel2str = f"(chainID A B C D E and {polarAtoms} and not (resid {residue} and chainID {chain}))"

                for _ in tqdm(u.trajectory):

                    resTupleSet = set()
                    for atom in u.select_atoms(f"{sel2str} and around {cutoff} {sel1str}"):

                        atomTuple = (atom.resid, atom.chainID, atom.name)
                        resTuple = (atom.resid, atom.chainID)

                        if atomTuple not in atomDict:
                            atomDict[atomTuple] = 1
                        else:
                            atomDict[atomTuple] += 1

                        resTupleSet.add(resTuple)

                    for element in resTupleSet:
                        if element not in residueDict:
                            residueDict[element] = 1
                        else:
                            residueDict[element] += 1

                    nFrames += 1

            # print(residueDict)  # debug
            # print(atomDict)     # debug

            # Make the ATOM occupancies fractional by dividing by the amount of hits
            # the residue received in total.
            for key1 in residueDict:
                for key2 in atomDict:
                    if key1[0] == key2[0] and key1[1] == key2[1]:
                        atomDict[key2] /= float(residueDict[key1])

            # Make the RESIDUE occupancies fractional by dividing the number of hits
            # by the total number of frames.
            for key in residueDict:
                residueDict[key] /= float(nFrames)

            # Sort the the residueDict by value
            residueDict = dict(sorted(residueDict.items(), key=lambda item: item[1], reverse=True))

            # print(residueDict)  # debug
            # print(atomDict)     # debug

            # Write to file.
            with open(f"newcontacts/{sim}_{residue}_{chain}.txt", 'w') as file:
                for key in residueDict:
                    file.write(f"{key[0]:<4d} {key[1]} {residueDict[key]:.4f}\n")
