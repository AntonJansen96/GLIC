#!/usr/bin/env python3

import MDAnalysis
import multiprocessing as mp
import numpy as np
from MDAnalysis.analysis.distances import distance_array
import matplotlib
import matplotlib.pyplot as plt

from science.parsing import loadCol
from science.utility import makeSuperDict
from science.utility import inputOptionHandler

# PARAMETERS

combs = [['E82', 'K38c']]

sims    = ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']
reps    = [1, 2, 3, 4]
chains  = ['A', 'B', 'C', 'D', 'E']
metrics = ['beta_expansion']

carboxylAtoms = 'name OE1 OE2 OD1 OD2 NE2 ND1'
polarAtoms    = 'name NZ NE2 OD2 OE1 OD1 NE OG1 OH ND2 OE2 OT1 NH2 O NH1 NE1 N ND1 OT2 OG'

# TASK FOR MULTIPROCESSING

def task(residue, target, sim, rep, chain1, chain2):

    def dist(v1, v2) -> float:
        x = v1[0] - v2[0]
        y = v1[1] - v2[1]
        z = v1[2] - v2[2]
        return (x * x + y * y + z * z)**0.5

    p1 = f'../../sims/{sim}/{rep:02d}/CA.pdb'
    p2 = f'../../sims/{sim}/{rep:02d}/MD_conv.xtc'
    u = MDAnalysis.Universe(p1, p2)

    superData = makeSuperDict([metrics, []])

    # * Take this out of the loop for speedup.

    if residue != 'all':
        sel1 = u.select_atoms(f'chainID {chain1} and resid {residue} and {carboxylAtoms}')
        sel2 = u.select_atoms(f'chainID {chain2} and resid {target} and {polarAtoms}')

    ECD_su = u.select_atoms(f"resid   0:193 and name CA and chainID {chain1}")
    ECD    = u.select_atoms( "resid   0:193 and name CA")
    TMD    = u.select_atoms( "resid 194:315 and name CA")
    TMD_su = u.select_atoms(f"resid 194:315 and name CA and chainID {chain1}")

    beta1 = u.select_atoms(f"chainID {chain1} and name CA and resid  30:34 ")
    beta2 = u.select_atoms(f"chainID {chain1} and name CA and resid 190:194")

    chainMap = {'A': 'E', 'B': 'A', 'C': 'B', 'D': 'C', 'E': 'D'}
    m2 = u.select_atoms(f"resid 241:245 and name CA and chainID {chain1}")
    m1 = u.select_atoms(f"resid 200:204 and name CA and chainID {chainMap[chain1]}")

    ca = u.select_atoms('resid 233 and name CA')
    ninePrime = u.select_atoms(f'resid 233 and chainID {chain1}')

    loopc = u.select_atoms(f'resid 174:182 and name CA and chainID {chain1}')
    xloop = u.select_atoms(f'resid 134:140 and name CA and chainID {chain1}')
    yloop = u.select_atoms(f'resid 145:150 and name CA and chainID {chainMap[chain1]}')

    # Start looping through the trajectory.
    for _ in u.trajectory:
        # If our selection is making a contact, continue with the analysis.
        if (residue == 'all') or (distance_array(sel1, sel2).min() < 4.0):

            if 'ecd_twist' in metrics:

                ECD_su_com = ECD_su.center_of_mass()
                ECD_com    = ECD.center_of_mass()
                TMD_com    = TMD.center_of_mass()
                TMD_su_com = TMD_su.center_of_mass()
                ecd_twist_coords   = np.array([ECD_su_com, ECD_com, TMD_com, TMD_su_com])
                ecd_twist_universe = MDAnalysis.Universe.empty(4, trajectory=True)
                ecd_twist_universe.atoms.positions = ecd_twist_coords
                superData['ecd_twist'].append(MDAnalysis.core.topologyobjects.Dihedral([0, 1, 2, 3], ecd_twist_universe).dihedral())

            if 'beta_expansion' in metrics:

                beta_com1 = beta1.center_of_mass()
                beta_com2 = beta2.center_of_mass()
                superData['beta_expansion'].append(dist(beta_com1, beta_com2))

            if 'm2_m1_dist' in metrics:

                m2_com = m2.center_of_mass()
                m1_com = m1.center_of_mass()
                superData['m2_m1_dist'].append(dist(m1_com, m2_com))

            if 'nine_prime_dist' in metrics:

                ca_com = ca.center_of_mass()
                min_dist = 10000000
                for pos in ninePrime.positions:
                    distance = dist(pos, ca_com)
                    if distance < min_dist:
                        min_dist = distance
                superData['nine_prime_dist'].append(min_dist)

            if 'loopc' in metrics:
                loopc_com = loopc.center_of_mass()
                xloop_com = xloop.center_of_mass()
                superData['loopc'].append(dist(loopc_com, xloop_com))

            if 'loopc2' in metrics:
                loopc_com = loopc.center_of_mass()
                yloop_com = yloop.center_of_mass()
                superData['loopc2'].append(dist(loopc_com, yloop_com))

    with open(f'{residue}_{target}_{sim}_{rep}_{chain1}_{chain2}.txt', 'w') as file:
        # Write header.
        for metric in metrics:
            file.write('{} '.format(metric))
        file.write('\n')
        # Write data.
        for kk in range(0, len(superData[metrics[0]])):
            for metric in metrics:
                file.write('{:.4f} '.format(superData[metric][kk]))
            file.write('\n')

# MAIN CODE

opMode = inputOptionHandler('Choose option', ['do nothing', 'run data', 'run plot'])

if opMode == 1:

    for comb in combs:

        # PREPROCESSING

        fullResidueName = comb[0]           # E35
        fullTargetName  = comb[1]           # T158c

        if fullResidueName == 'all':
            residue = fullResidueName
            target  = fullTargetName
        else:
            residue = int(fullResidueName[1:])  # 35
            target  = fullTargetName[1:]        # 158c

        if target[-1] == 'p':
            chains2 = ['B', 'C', 'D', 'E', 'A']
            target = target[:-1]  # Remove trailing letter (158c -> 158).
        elif target[-1] == 'c':
            chains2 = ['E', 'A', 'B', 'C', 'D']
            target = target[:-1]  # Remove trailing letter (158c -> 158).
        else:
            chains2 = ['A', 'B', 'C', 'D', 'E']

        # GATHER ITERABLES

        items = []
        for sim in sims:
            for rep in reps:
                for ii in range(0, len(chains)):
                    items.append((residue, target, sim, rep, chains[ii], chains2[ii]))

        # RUN MULTIPROCESSING

        pool = mp.Pool(processes=mp.cpu_count())
        pool.starmap(task, items, chunksize=1)

if opMode == 2:

    def stderr(array):
        return np.std(array) / np.sqrt(len(array))

    matplotlib.rcParams.update({'font.size': 20})

    for comb in combs:

        # PREPROCESSING

        fullResidueName = comb[0]           # E35
        fullTargetName  = comb[1]           # T158c
        residue = int(fullResidueName[1:])  # 35
        target  = fullTargetName[1:]        # 158c

        if target[-1] == 'p':
            chains2 = ['B', 'C', 'D', 'E', 'A']
            target = target[:-1]  # Remove trailing letter (158c -> 158).
        elif target[-1] == 'c':
            chains2 = ['E', 'A', 'B', 'C', 'D']
            target = target[:-1]  # Remove trailing letter (158c -> 158).
        else:
            chains2 = ['A', 'B', 'C', 'D', 'E']

        # GATHER DATA / full trajectories from 'bloom' directory.

        fullData   = makeSuperDict([sims, reps, chains, metrics, []])
        fullResult = makeSuperDict([sims, metrics, []])

        for sim in sims:
            for rep in reps:
                for chain in chains:
                    fname = f'all_all_{sim}_{rep}_{chain}_{chain}.txt'
                    for idx in range(0, len(metrics)):
                        fullData[sim][rep][chain][metrics[idx]] = loadCol(fname, idx + 1, header=0)

        for sim in sims:
            for metric in metrics:
                for rep in reps:
                    for chain in chains:
                        fullResult[sim][metric].append(np.mean(fullData[sim][rep][chain][metric]))

        # GATHER DATA / filtered data from the contact.

        filteredData   = makeSuperDict([sims, reps, chains, metrics, []])
        filteredResult = makeSuperDict([sims, metrics, []])

        for sim in sims:
            for rep in reps:
                for ii in range(0, len(chains)):
                    fname = f'{residue}_{target}_{sim}_{rep}_{chains[ii]}_{chains2[ii]}.txt'
                    for idx in range(0, len(metrics)):
                            filteredData[sim][rep][chains[ii]][metrics[idx]] = loadCol(fname, idx + 1, header=0)

        for sim in sims:
            for metric in metrics:
                for rep in reps:
                    for chain in chains:
                        # If a file is empty this will results in a np.mean being
                        # taken over an empty list, which will results in a np.NAN
                        # being returned. We need to remove these before we proceed.
                        array = filteredData[sim][rep][chain][metric]
                        if len(array) > 0:
                            filteredResult[sim][metric].append(np.mean(array))

        # MAKE THE BARPLOTS

        nameList = ['All', f'{fullResidueName}-{fullTargetName}']
        width    = 0.2
        x        = np.arange(len(nameList))

        for metric in metrics:

            meanList = {'4HFI_4': [], '4HFI_7': [], '6ZGD_4': [], '6ZGD_7': []}
            serrList = {'4HFI_4': [], '4HFI_7': [], '6ZGD_4': [], '6ZGD_7': []}
            for sim in sims:
                meanList[sim].append(np.mean(fullResult[sim][metric]))
                serrList[sim].append(stderr(fullResult[sim][metric]))
                meanList[sim].append(np.mean(filteredResult[sim][metric]))
                serrList[sim].append(stderr(filteredResult[sim][metric]))

            fig = plt.figure(figsize=(8, 6))
            ax  = fig.add_subplot()

            # 6ZGD_7
            mean4 = meanList['6ZGD_7']
            serr4 = serrList['6ZGD_7']
            ax.bar(     x - width * 1.5, mean4, width, color='#8856a7', label='closed, pH 7')
            ax.errorbar(x - width * 1.5, mean4, serr4, color='#8856a7', fmt='none', capsize=6, linewidth=2)

            # 6ZGD_4
            mean3 = meanList['6ZGD_4']
            serr3 = serrList['6ZGD_4']
            ax.bar(     x - width / 2.0, mean3, width, color='#9ebcda', label='closed, pH 4')
            ax.errorbar(x - width / 2.0, mean3, serr3, color='#9ebcda', fmt='none', capsize=6, linewidth=2)

            # 4HFI_7
            mean2 = meanList['4HFI_7']
            serr2 = serrList['4HFI_7']
            ax.bar(     x + width / 2.0, mean2, width, color='#8856a7', label='open, pH 7', edgecolor='w', lw=1, hatch='//')
            ax.errorbar(x + width / 2.0, mean2, serr2, color='#8856a7', fmt='none', capsize=6, linewidth=2)

            # 4HFI_4
            mean1 = meanList['4HFI_4']
            serr1 = serrList['4HFI_4']
            ax.bar(     x + width * 1.5, mean1, width, color='#9ebcda', label='open, pH 4', edgecolor='w', lw=1, hatch='\\\\')
            ax.errorbar(x + width * 1.5, mean1, serr1, color='#9ebcda', fmt='none', capsize=6, linewidth=2)

            ax.set_xticks(x, nameList)
            ax.legend(loc=1, prop={'size': 12})

            if metric == 'ecd_twist':
                plt.ylabel('ECD Twist (degrees)')
                plt.ylim(-12.5, -19)

            if metric == 'beta_expansion':
                plt.ylabel('Beta Expansion (Å)')
                plt.ylim(14, 17)

            if metric == 'm2_m1_dist':
                plt.ylabel('M2-M1(-) Distance (Å)')
                plt.ylim(11, 20)

            if metric == 'nine_prime_dist':
                plt.ylabel('9\' Radius (Å)')
                plt.ylim(2.5, 4)

            if metric == 'loopc':
                plt.ylabel('loopC-loopB Distance (Å)')
                plt.ylim(9, 16)

            if metric == 'loopc2':
                plt.ylabel('loopC-comp Distance (Å)')
                plt.ylim(13, 19)

            plt.legend()
            outname = f'{residue}_{target}_{metric}.png'
            plt.tight_layout()
            plt.savefig(outname)
            plt.clf()
            plt.close()

        # MAKE THE HISTOGRAMS

        for metric in metrics:

            # Define metric-specific stuff

            if metric == 'ecd_twist':
                histRange = (-25, -5)
                xlabel = 'ECD Twist (degrees)'

            if metric == 'beta_expansion':
                histRange = (12, 18)
                xlabel = 'Beta Expansion (Å)'

            if metric == 'm2_m1_dist':
                histRange = (9, 21)
                xlabel = 'M2-M1(-) Distance (Å)'

            if metric == 'nine_prime_dist':
                histRange = (2, 5.5)
                xlabel = '9\' Radius (Å)'

            if metric == 'loopc':
                histRange = (8, 18)
                xlabel = 'loopC-loopB Distance (Å)'

            if metric == 'loopc2':
                histRange = (12, 20)
                xlabel = 'loopC-comp Distance (Å)'

            # Start building the figure

            fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(8, 8))

            for idx in [0, 1]:

                maxval = []
                subplt = axs[idx]

                for sim in sims:

                    # Each figure contains four curves (one for each sim).

                    histList = []

                    for rep in reps:
                        for chain in chains:

                            # Get histogram values and corresponding bin_edges
                            if idx == 0:
                                x = fullData[sim][rep][chain][metric]
                            if idx == 1:
                                x = filteredData[sim][rep][chain][metric]

                            hist, bin_edges = np.histogram(x, density=True, bins=30, range=histRange)

                            # Another bug fix to prevent nan from np.mean([])
                            if 'nan' in [str(val) for val in hist]:
                                continue

                            histList.append(hist)

                    # Compute mean and standard error.
                    meanList = len(hist) * [0]
                    errorList = len(hist) * [0]

                    for kk in range(0, len(hist)):

                        # Create list of 20 values
                        temp = [0] * len(histList)  # 4*5=20
                        for ll in range(0, len(histList)):  # 4*5=20
                            temp[ll] = histList[ll][kk]

                        meanList[kk]  = np.mean(temp)
                        errorList[kk] = stderr(temp)

                    # Plot mean and shaded region (standard error).
                    A = []
                    B = []
                    for kk in range(0, len(meanList)):
                        A.append(meanList[kk] + errorList[kk])
                        B.append(meanList[kk] - errorList[kk])

                    # The values in these dictionaries come from bloom/static.txt, which
                    # in tern come from 4HFI_clean_crystal.pdb and 6ZGD_clean_cryoEM.pdb.
                    # Margins in matplotlib are 1.05 * largest value = max(maxval).
                    # This is to make sure that the vlines reach top of figure without rescaling.
                    static4HFI = {'ecd_twist': -11.631, 'beta_expansion': 13.767, 'm2_m1_dist': 14.094, 'nine_prime_dist': 5.075, 'loopc': 0, 'loopc2': 0}
                    static6ZGD = {'ecd_twist': -16.832, 'beta_expansion': 15.430, 'm2_m1_dist': 17.270, 'nine_prime_dist': 3.497, 'loopc': 0, 'loopc2': 0}

                    if sim == '4HFI_4':
                        subplt.plot(bin_edges[1:], meanList, linewidth=2,   color='#9ebcda', label='open, pH 4', linestyle='--')
                        subplt.fill_between(bin_edges[1:], A, B, alpha=0.5, color='#9ebcda', hatch='\\\\', edgecolor='w', lw=1)

                    if sim == '4HFI_7':
                        subplt.vlines(x=static4HFI[metric], ymin=0, ymax=10, color='black', label='open, crystal', linestyle='--')
                        subplt.plot(bin_edges[1:], meanList, linewidth=2,   color='#8856a7', label='open, pH 7', linestyle='--')
                        subplt.fill_between(bin_edges[1:], A, B, alpha=0.5, color='#8856a7', hatch='//',   edgecolor='w', lw=1)

                    if sim == '6ZGD_4':
                        subplt.plot(bin_edges[1:], meanList, linewidth=2,   color='#9ebcda', label='closed, pH 4')
                        subplt.fill_between(bin_edges[1:], A, B, alpha=0.5, color='#9ebcda')

                    if sim == '6ZGD_7':
                        subplt.vlines(x=static6ZGD[metric], ymin=0, ymax=10, color='black', label='closed, cryoEM')
                        subplt.plot(bin_edges[1:], meanList, linewidth=2,   color='#8856a7', label='closed, pH 7')
                        subplt.fill_between(bin_edges[1:], A, B, alpha=0.5, color='#8856a7')

                    maxval.append(max(A))

                subplt.set_ylim(0, max(maxval) * 1.05)
                subplt.set_xlim(histRange)

                if idx == 0:
                    subplt.set_ylabel('All')
                    subplt.legend(loc='best', prop={'size': 14})

                if idx == 1:
                    subplt.set_ylabel(f'{fullResidueName}-{fullTargetName}')
                    subplt.set_xlabel(xlabel)

            outname = f'{residue}_{target}_{metric}_hist.png'
            plt.tight_layout()
            plt.savefig(outname)
            plt.clf()
            plt.close()
