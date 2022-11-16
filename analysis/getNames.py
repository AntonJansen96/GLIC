#!/bin/python3

import MDAnalysis

u       = MDAnalysis.Universe('../sims/4HFI_4/01/CA.pdb')
protein = u.select_atoms('chainID A to E')

uniqueNames = set()
uniqueResnames = set()

for resname in list(protein.residues.resnames):
    uniqueResnames.add(resname)

for atomname in list(protein.atoms.names):
    if ('O' in atomname) or ('N' in atomname):
        uniqueNames.add(atomname)

print(uniqueResnames)
print(uniqueNames)
