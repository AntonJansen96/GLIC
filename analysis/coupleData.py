#!/usr/bin/env python3

import MDAnalysis
import multiprocessing as mp
import numpy as np
from MDAnalysis.analysis.distances import distance_array
from science.utility import makeSuperDict

# PARAMETERS

alla  = [['all', 'all']]
regul = [['E26', 'V79p'], ['E26', 'N80p'], ['E26', 'V81p'], ['E82', 'T36'], ['E82', 'K38c'], ['E35', 'L114'], ['E243', 'K248'], ['E243', 'N200c']]
loopC = [['E177', 'K148c'], ['D178', 'K148c'], ['E181', 'R133'], ['D185', 'I128']]
combs = alla + regul + loopC

sims    = ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']
reps    = [1, 2, 3, 4]
chains  = ['A', 'B', 'C', 'D', 'E']
metrics = ['ecd_twist', 'beta_expansion', 'm2_m1_dist', 'nine_prime_dist', 'loopc', 'loopc2']

carboxylAtoms = 'name OE1 OE2 OD1 OD2 NE2 ND1'
polarAtoms    = 'name NZ NE2 OD2 OE1 OD1 NE OG1 OH ND2 OE2 OT1 NH2 O NH1 NE1 N ND1 OT2 OG'

# TASK FOR MULTIPROCESSING

def task(residue, target, sim, rep, chain1, chain2):

    def dist(v1, v2) -> float:
        x = v1[0] - v2[0]
        y = v1[1] - v2[1]
        z = v1[2] - v2[2]
        return (x * x + y * y + z * z)**0.5

    p1 = f'../sims/{sim}/{rep:02d}/CA.pdb'
    p2 = f'../sims/{sim}/{rep:02d}/MD_conv.xtc'
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

    with open(f'couple/{residue}_{target}_{sim}_{rep}_{chain1}_{chain2}.txt', 'w') as file:
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
