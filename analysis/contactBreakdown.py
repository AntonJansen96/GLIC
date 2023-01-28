#!/usr/bin/env python3

import MDAnalysis
import multiprocessing as mp

from science.utility import createIndexFile
from science.utility import gromacs
from science.utility import makeSuperDict
from science.parsing import loadxvg

if __name__ == "__main__":

    for comb in [['E243', 'K248']]:

        # Some name variables
        fullResidueName = comb[0]           # E35
        fullTargetName  = comb[1]           # T158c
        residue = int(fullResidueName[1:])  # 35
        target  = fullTargetName[1:]        # 158c

        # Select the correct chain, depending on principal, complementary, etc.
        chain1 = ['A', 'B', 'C', 'D', 'E']
        if target[-1] == 'p':
            chain2 = ['B', 'C', 'D', 'E', 'A']
            target = target[:-1]  # Remove the trailing letter (158c -> 158).
        elif target[-1] == 'c':
            chain2 = ['E', 'A', 'B', 'C', 'D']
            target = target[:-1]  # Remove the trailing letter (158c -> 158).
        else:
            chain2 = ['A', 'B', 'C', 'D', 'E']

        # The sims and reps we loop over.
        sims = ['6ZGD_7', '6ZGD_4', '4HFI_7', '4HFI_4']
        reps = [1, 2, 3, 4]

        # Get the NAMES of the atoms in the target (T158c) that we need to consider.
        # We basically just want the NAMES of the polar atoms in the target (T158c).
        # For this we can simply look hardcoded in chain A as all five chains are the same.
        # We can also simply use replica one of 4HFI_4 for ALL as the residues are all the same.
        temp  = f'resid {target} and chainID A and name NZ NE2 OD2 OE1 OD1 NE OG1 OH ND2 OE2 OT1 NH2 O NH1 NE1 N ND1 OT2 OG'
        u     = MDAnalysis.Universe('../sims/4HFI_4/01/CA.pdb')
        names = list(u.select_atoms(temp).names)
        # print(names)  # debug

        #! Task for multithreading.
        def task(p1, indexFileName, sel1, sel2, reps, sim, fullResidueName, fullTargetName, chain1IDX, chain2IDX, name):
            createIndexFile(p1, indexFileName, [sel1, sel2])

            # Loop over each replica and run GROMACS mindist
            for rep in reps:
                p2 = f'../sims/{sim}/{rep:02d}/MD_conv.xtc'
                xvgFileName = f'backbone/{sim}_{rep}_{fullResidueName}_{fullTargetName}_{chain1IDX}{chain2IDX}_{name}.xvg'
                gromacs(f'mindist -s {p1} -f {p2} -n {indexFileName} -od {xvgFileName}', stdin=[0, 1])

        #! Gather iterables.
        items = []

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

                    # THIS BLOCK OF CODE WE COMMENT OUT AS WE WANT TO MULTITHREAD IT ####
                    # createIndexFile(p1, indexFileName, [sel1, sel2])
                    # # Loop over each replica and run GROMACS mindist
                    # for rep in reps:
                    #     p2 = f'../sims/{sim}/{rep:02d}/MD_conv.xtc'
                    #     xvgFileName = f'backbone/{sim}_{rep}_{fullResidueName}_{fullTargetName}_{chain1[idx]}{chain2[idx]}_{name}.xvg'
                    #     gromacs(f'mindist -s {p1} -f {p2} -n {indexFileName} -od {xvgFileName}', stdin=[0, 1])
                    # END OF BLOCK OF CODE WE COMMENT OUT AS WE WANT TO MULTITHREAD IT ####

                    items.append((p1, indexFileName, sel1, sel2, reps, sim, fullResidueName, fullTargetName, chain1[idx], chain2[idx], name))

        #! Run multithreaded.
        pool = mp.Pool(processes=mp.cpu_count())
        pool.starmap(task, items, chunksize=1)

        #? Gather data and process.
        superData = makeSuperDict([sims, names, 0])
        superDataPercent = makeSuperDict([sims, names, 0])

        for sim in sims:

            totalCount = 0  # total #frames with a contact (for all atoms) in all reps and chains.

            for name in names:

                count = 0   # total #frames with a contact for a specific atom in all reps and chains.

                for rep in reps:
                    for idx in range(0, len(chain1)):

                        x = loadxvg(f'backbone/{sim}_{rep}_{fullResidueName}_{fullTargetName}_{chain1[idx]}{chain2[idx]}_{name}.xvg')[1]

                        for dist in x:
                            if dist < 0.40:  # cutoff is 0.4nm.
                                count += 1
                                totalCount += 1

                superData[sim][name] = count

            for name in names:
                if totalCount == 0:
                    superDataPercent[sim][name]
                else:
                    superDataPercent[sim][name] = round(superData[sim][name] / float(totalCount), 2)

        # print(superData)         # debug
        # print(superDataPercent)  # debug

        # Neatly write the result to an output file.
        with open('contactBreakdown.txt', 'a+') as file:
            file.write(f'{fullResidueName}-{fullTargetName}:\n')
            for sim in sims:

                totalCount = 0
                for name in names:
                    totalCount += superData[sim][name]
                totalPercent = 100 * (totalCount / (1000.0 * 5 * len(reps)))

                for name in names:
                    file.write(f'{sim} {name:4s}: {superDataPercent[sim][name]:.2f} {superData[sim][name]:7d}   (total = {totalCount:5d}, {totalPercent:.1f}% occupancy)\n')
            file.write('\n')
