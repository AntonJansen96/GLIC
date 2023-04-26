#!/usr/bin/env python3

import MDAnalysis
from tqdm import tqdm
from science.utility import triplet2letter

# PARAMETERS #############################################################################

sims = ['6ZGD_7']
reps = [1]
residues = [35]
chains = ['A']

carboxylAtoms = 'name OE1 OE2 OD1 OD2 NE2 ND1'
polarAtoms = 'name NZ NE2 OD2 OE1 OD1 NE OG1 OH ND2 OE2 OT1 NH2 O NH1 NE1 N ND1 OT2 OG'

cutoff = 4.0  # Angstroms.

##########################################################################################

for residue in residues:
    for sim in sims:

        nFrames = 0
        atomDict = {}
        residueDict = {}

        for chain in chains:
            for rep in reps:

                path1 = f"../sims/{sim}/{rep:02d}/CA.pdb"
                path2 = f"../sims/{sim}/{rep:02d}/MD_conv.xtc"
                u = MDAnalysis.Universe(path1, path2)

                target = f"(resid {residue} and chainID {chain} and {carboxylAtoms})"
                consider = f"(chainID A B C D E and {polarAtoms} and not (resid {residue} and chainID {chain}))"

                for _ in tqdm(u.trajectory):

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

        print(residueDict)  # debug
        print(atomDict)     # debug

        # Make the ATOM occupancies fractional by dividing by the amount of hits
        # the residue received in total.
        for key1 in residueDict:
            for key2 in atomDict:
                if key1 == key2[0]:
                    atomDict[key2] /= float(residueDict[key1])

        # Make the RESIDUE occupancies fractional by dividing the number of hits
        # by the total number of frames.
        for key in residueDict:
            residueDict[key] /= float(nFrames)

        # Sort the the residueDict by value
        residueDict = dict(sorted(residueDict.items(), key=lambda item: item[1], reverse=True))

        print(residueDict)  # debug
        print(atomDict)     # debug

        # Write to file.
        with open(f"newcontacts/{sim}_{residue}.txt", 'w') as file:
            for key in residueDict:
                file.write(f"{key:<5s} {residueDict[key]:.3f}\n")

            file.write('\n')
            for key1 in residueDict:
                for key2 in atomDict:
                    if key1 == key2[0]:
                        file.write(f"{key1:<5s} {key2[1]:<4s} {atomDict[key2]:.3f}\n")
