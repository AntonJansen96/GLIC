#!/usr/bin/env python3

import MDAnalysis
from MDAnalysis.analysis.distances import distance_array
import multiprocessing as mp
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from science.utility import makeSuperDict
from science.parsing import loadCol

# Set global font size for figures.
matplotlib.rcParams.update({'font.size': 20})

def stderr(array):
    return np.std(array) / np.sqrt(len(array))

# PARAMETERS.

sims      = ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']
reps      = [1, 2, 3, 4]
chains1   = ['A', 'B', 'C', 'D', 'E']

carboxylAtoms = 'name OE1 OE2 OD1 OD2 NE2 ND1'
polarAtoms    = 'name NZ NE2 OD2 OE1 OD1 NE OG1 OH ND2 OE2 OT1 NH2 O NH1 NE1 N ND1 OT2 OG'

# for comb in [['E26', 'V79p'], ['E26', 'N80p'], ['E26', 'V81p'], ['E82', 'T36'], ['E82', 'K38c'], ['E35', 'L114'], ['E243', 'K248'], ['E243', 'N200c']]:
for comb in [['E26', 'V79p']]:

    # PROCESS NAMES.
    fullResidueName = comb[0]           # E35
    fullTargetName  = comb[1]           # T158c
    residue = int(fullResidueName[1:])  # 35
    target  = fullTargetName[1:]        # 158c

    # SELECT CORRECT CHAIN.
    if target[-1] == 'p':
        chains2 = ['B', 'C', 'D', 'E', 'A']
        target = target[:-1]  # Remove the trailing letter (158c -> 158).
    elif target[-1] == 'c':
        chains2 = ['E', 'A', 'B', 'C', 'D']
        target = target[:-1]  # Remove the trailing letter (158c -> 158).
    else:
        chains2 = ['A', 'B', 'C', 'D', 'E']

    # MULTIPROCESSING TASK.
    def task(sim, rep, residue, target, chains1, chains2):
        fname = f'couple/{residue}_{target}_{sim}_{rep}_{chains1[0]}_{chains2[0]}.txt'
        if not os.path.isfile(fname):

            p1 = f'../sims/{sim}/{rep:02d}/CA.pdb'
            p2 = f'../sims/{sim}/{rep:02d}/MD_conv.xtc'
            u  = MDAnalysis.Universe(p1, p2)

            for idx in range(0, len(chains1)):

                sel1 = u.select_atoms(f'chainID {chains1[idx]} and resid {residue} and {carboxylAtoms}')
                sel2 = u.select_atoms(f'chainID {chains2[idx]} and resid {target}  and {polarAtoms}')

                fname = f'couple/{residue}_{target}_{sim}_{rep}_{chains1[idx]}_{chains2[idx]}.txt'
                with open(fname, 'w+') as file:
                    for ts in u.trajectory:
                        if distance_array(sel1, sel2).min() < 4.0:
                            file.write(f'{ts.frame}\n')

    # PREPARE ITERABLES.
    items = []
    for sim in sims:
        for rep in reps:
            items.append((sim, rep, residue, target, chains1, chains2))

    # RUN MULTITHREADED.
    pool = mp.Pool(processes=mp.cpu_count())
    pool.starmap(task, items, chunksize=1)

    #! BLOOM DATA GENERATION PART ##############################################

    # TASK
    def doBlooming(sim, rep, residue, target, chain1, chain2):

        def dist(c1, c2) -> float:
            x = c1[0] - c2[0]
            y = c1[1] - c2[1]
            z = c1[2] - c2[2]
            return (x * x + y * y + z * z)**0.5

        fname = f'couple/{residue}_{target}_{sim}_{rep}_{chain1}_{chain2}.dat'
        if not os.path.exists(fname):

            p1 = f'../sims/{sim}/{rep:02d}/CA.pdb'
            p2 = f'../sims/{sim}/{rep:02d}/MD_conv.xtc'
            u  = MDAnalysis.Universe(p1, p2)

            for ii in range(0, len(chains1)):

                metrics = ['ecd_twist', 'beta_expansion', 'm2_m1_dist', 'nine_prime_dist']
                superData = makeSuperDict([metrics, []])

                fname = f'couple/{residue}_{target}_{sim}_{rep}_{chain1}_{chain2}.txt'
                contactFrameIndices = [int(val) for val in loadCol(fname)]

                for jj in contactFrameIndices:
                    u.trajectory[jj]

                    # *** ecd_twist ***

                    ECD_su_com = u.select_atoms(f"protein and resid   0:193 and name CA and chainID {chain1}").center_of_mass()
                    ECD_com    = u.select_atoms( "protein and resid   0:193 and name CA").center_of_mass()
                    TMD_com    = u.select_atoms( "protein and resid 194:315 and name CA").center_of_mass()
                    TMD_su_com = u.select_atoms(f"protein and resid 194:315 and name CA and chainID {chain1}").center_of_mass()
                    ecd_twist_coords   = np.array([ECD_su_com, ECD_com, TMD_com, TMD_su_com])
                    ecd_twist_universe = MDAnalysis.Universe.empty(4, trajectory=True)
                    ecd_twist_universe.atoms.positions = ecd_twist_coords
                    superData['ecd_twist'].append(MDAnalysis.core.topologyobjects.Dihedral([0, 1, 2, 3], ecd_twist_universe).dihedral())

                    # *** beta_expansion ***

                    beta_com1 = u.select_atoms(f"protein and chainID {chain1} and name CA and resid  30:34 ").center_of_mass()
                    beta_com2 = u.select_atoms(f"protein and chainID {chain1} and name CA and resid 190:194").center_of_mass()
                    superData['beta_expansion'].append(dist(beta_com1, beta_com2))

                    # *** m2_m1_dist ***

                    m2_com = u.select_atoms(f"protein and resid 241:245 and name CA and chainID {chain1}").center_of_mass()
                    m1_com = u.select_atoms(f"protein and resid 200:204 and name CA and chainID {chain2}").center_of_mass()
                    superData['m2_m1_dist'].append(dist(m1_com, m2_com))

                    # *** nine_prime_dist ***

                    min_dist = 10000000
                    resid1 = u.select_atoms(f'protein and resid 233 and chainID {chain1}')
                    ca_com = u.select_atoms( 'protein and resid 233 and name CA').center_of_mass()
                    for pos in resid1.positions:
                        distance = dist(pos, ca_com)
                        if distance < min_dist:
                            min_dist = distance
                    superData['nine_prime_dist'].append(min_dist)

                    # *** Anton: come up with your own metric for loop-C here, e.g. RMSD ***

                with open(f'couple/{residue}_{target}_{sim}_{rep}_{chain1}_{chain2}.dat', 'w') as file:
                    # Write header.
                    for metric in metrics:
                        file.write('{} '.format(metric))
                    file.write('\n')
                    # Write data.
                    for kk in range(0, len(superData['nine_prime_dist'])):
                        for metric in metrics:
                            file.write('{:.4f} '.format(superData[metric][kk]))
                        file.write('\n')

    # GATHER ITERABLES
    items = []
    for sim in sims:
        for rep in reps:
            for idx in range(0, len(chains1)):
                items.append((sim, rep, residue, target, chains1[idx], chains2[idx]))

    # RUN MULTITHREADED
    pool = mp.Pool(processes=mp.cpu_count())
    pool.starmap(doBlooming, items, chunksize=1)

    #! VISUALIZATION PART ######################################################

    #? Load+process the data for the full trajectories from bloom.

    metrics    = ['ecd_twist', 'ecd_spread', 'ecd_upper_spread', 'beta_expansion', 'm2_m1_dist', 'm2_radius', 'm1_kink', 'm1_kink_alt', 'nine_prime_dist', 'nine_prime_pore', 'minus_two_prime_dist', 'minus_two_prime_pore', 'c_loop']
    fullData   = makeSuperDict([sims, reps, chains1, metrics, []])
    fullResult = makeSuperDict([sims, metrics, []])

    # Load the data.
    for sim in sims:
        for rep in reps:
            for chain in chains1:
                fname = f'bloom/{sim}_{rep}_{chain}.txt'
                for idx in range(0, len(metrics)):
                    fullData[sim][rep][chain][metrics[idx]] = loadCol(fname, idx + 1, header=0)

    # Get the results.
    for sim in sims:
        for metric in metrics:
            for rep in reps:
                for chain in chains1:
                    fullResult[sim][metric].append(np.mean(fullData[sim][rep][chain][metric]))

    #? Load+process the data filtered for the contact.

    metrics        = ['ecd_twist', 'beta_expansion', 'm2_m1_dist', 'nine_prime_dist']
    filteredData   = makeSuperDict([sims, reps, chains1, metrics, []])
    filteredResult = makeSuperDict([sims, metrics, []])

    # Load the data.
    for sim in sims:
        for rep in reps:
            for ii in range(0, len(chains1)):
                fname = f'couple/{residue}_{target}_{sim}_{rep}_{chains1[ii]}_{chains2[ii]}.dat'
                for idx in range(0, len(metrics)):
                        filteredData[sim][rep][chains1[ii]][metrics[idx]] = loadCol(fname, idx + 1, header=0)

    # Get the results.
    for sim in sims:
        for metric in metrics:
            for rep in reps:
                for chain in chains1:
                    # If a file is empty this will results in a np.mean being
                    # taken over an empty list, which will results in a np.NAN
                    # being returned. We need to remove these before we proceed.
                    array = filteredData[sim][rep][chain][metric]
                    if len(array) > 0:
                        filteredResult[sim][metric].append(np.mean(array))

    #? Make BARPLOTS, comparing the two sets.

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

        ['ecd_twist', 'beta_expansion', 'm2_m1_dist', 'nine_prime_dist']

        if metric == 'ecd_twist':
            plt.ylabel('ECD Twist (degrees)')
            plt.ylim(-12.5, -19)

        if metric == 'beta_expansion':
            plt.ylabel('Beta Expansion (Å)')
            plt.ylim(14, 17)

        if metric == 'm2_m1_dist':
            plt.ylabel('M2-M1(-) Distance (Å)')
            plt.ylim(10, 18)

        if metric == 'nine_prime_dist':
            plt.ylabel('9\' Radius (Å)')
            plt.ylim(2.5, 4)

        plt.legend()
        outname = f'couple/{residue}_{target}_{metric}.png'
        plt.savefig(outname)
        os.system(f'convert {outname} -trim {outname}')
        plt.clf()
        plt.close()

    #? Make HISTOGRAMS, comparing the two sets (All vs contact).

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

        # Start building the figure

        fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(8, 8))

        for idx in [0, 1]:

            maxval = []
            subplt = axs[idx]

            for sim in sims:

                # Each figure contains four curves (one for each sim).

                histList = []

                for rep in reps:
                    for chain in chains1:

                        # Get histogram values and corresponding bin_edges
                        if idx == 0:
                            x = fullData[sim][rep][chain][metric]
                        if idx == 1:
                            x = filteredData[sim][rep][chain][metric]

                        hist, bin_edges = np.histogram(x, density=True, bins=20, range=histRange)

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
                static4HFI = {'ecd_twist': -11.631, 'ecd_spread': 23.796, 'ecd_upper_spread': 25.487, 'beta_expansion': 13.767, 'm2_m1_dist': 14.094, 'm2_radius': 21.583, 'm1_kink': 156.102, 'm1_kink_alt': 164.330, 'nine_prime_dist': 5.075, 'nine_prime_pore': 7.626, 'minus_two_prime_dist': 11.420, 'minus_two_prime_pore': 6.004, 'c_loop': 13.966}
                static6ZGD = {'ecd_twist': -16.832, 'ecd_spread': 23.548, 'ecd_upper_spread': 24.948, 'beta_expansion': 15.430, 'm2_m1_dist': 17.270, 'm2_radius': 19.425, 'm1_kink': 157.222, 'm1_kink_alt': 160.535, 'nine_prime_dist': 3.497, 'nine_prime_pore': 5.973, 'minus_two_prime_dist': 11.200, 'minus_two_prime_pore': 5.888, 'c_loop': 14.378}

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

        outname = f'couple/{residue}_{target}_{metric}_hist.png'
        plt.tight_layout()
        plt.savefig(outname)
        os.system(f'convert {outname} -trim {outname}')
        plt.clf()
        plt.close()
