#!/usr/bin/env python3

import sys
import MDAnalysis

from science.utility import createIndexFile
from science.utility import gromacs
from science.utility import makeSuperDict
from science.parsing import loadxvg


#? Handle CMDline input
assert len(sys.argv) == 3, "Specify residue and target, e.g. E35, T158c"
fullResidueName = sys.argv[1]
fullTargetName  = sys.argv[2]
residue         = int(sys.argv[1][1:])
target          = sys.argv[2][1:]
# print(fullResidueName, fullTargetName, residue, target)  # debug

#? Select chains depending on principal, complementary or same.
chain1 = ['A', 'B', 'C', 'D', 'E']
if target[-1] == 'p':
    chain2 = ['B', 'C', 'D', 'E', 'A']
elif target[-1] == 'c':
    chain2 = ['E', 'A', 'B', 'C', 'D']
else:
    chain2 = ['A', 'B', 'C', 'D', 'E']
target = target[:-1]
# print(chain1, chain2)  # debug

sims = ['4HFI_4']
reps = [1]

# Get names of the atoms in the target (T15c) that we need to consider.
temp  = f'resid {target} and chainID A and name NZ NE2 OD2 OE1 OD1 NE OG1 OH ND2 OE2 OT1 NH2 O NH1 NE1 N ND1 OT2 OG'
u     = MDAnalysis.Universe('../sims/4HFI_4/01/CA.pdb')
names = list(u.select_atoms(temp).names)
# print(names)  # debug

for sim in sims:

    for idx in range(0, len(chain1)):

        # Get MDAnalysis-style selection of residue (E35) carboxyl-oxygens in a specific chain.
        sel1 = f'resid {residue} and chainID {chain1[idx]} and name OE1 OE2 OD1 OD2 NE2 ND1'

        for name in names:
            sel2 = f'resid {target} and chainID {chain2[idx]} and name {name}'
            # print(f'{sel1} WITH {sel2}')  # debug

            # Make the index file (we only need one per four replicas)
            p1 = f'../sims/{sim}/01/CA.pdb'
            indexFileName = f'backbone/{sim}_{fullResidueName}_{fullTargetName}_{chain1[idx]}{chain2[idx]}_{name}.ndx'
            createIndexFile(p1, indexFileName, [sel1, sel2])

            # Loop over each replica and run GROMACS mindist
            for rep in reps:
                p2 = f'../sims/{sim}/{rep:02d}/MD_conv.xtc'
                xvgFileName = f'backbone/{sim}_{rep}_{fullResidueName}_{fullTargetName}_{chain1[idx]}{chain2[idx]}_{name}.xvg'
                gromacs(f'mindist -s {p1} -f {p2} -n {indexFileName} -od {xvgFileName}', stdin=[0, 1])

#? GATHER AND PROCESS DATA (this is fast enough)

superData = makeSuperDict([sims, names, 0])

for sim in sims:

    totalCount = 0  # total #frames with a contact (for all atoms) in all reps and chains

    for name in names:

        count = 0   # total #frames with a contact for a specific atom in all reps and chains

        for rep in reps:
            for idx in range(0, len(chain1)):

                x = loadxvg(f'backbone/{sim}_{rep}_{fullResidueName}_{fullTargetName}_{chain1[idx]}{chain2[idx]}_{name}.xvg')[1]

                for dist in x:
                    if dist < 0.4:
                        count += 1
                        totalCount += 1

        superData[sim][name] = count

    print(superData)

    for name in names:
        superData[sim][name] /= float(totalCount)

    print(superData)
