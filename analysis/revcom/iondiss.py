#!/usr/bin/env python3

import os
import numpy as np
import MDAnalysis
import multiprocessing as mp
from MDAnalysis.analysis.distances import distance_array

from science.cphmd import getLambdaFileIndices
from science.parsing import loadxvg
from science.parsing import loadCol

def stderr(array):
    return np.std(array) / np.sqrt(len(array))

# PARAMETERS ###################################################################

residues = [35]
# sims = ['4HFI_7']
sims = ['6ZGD_7', '6ZGD_4', '4HFI_7', '4HFI_4']
reps = [1, 2, 3, 4]
chains = ['A', 'B', 'C', 'D', 'E']
cutoff = 3.0  # nm

carboxylAtoms = "OE1 OE2 OD1 OD2 NE2 ND1"

################################################################################

# Task for multithreading.
def task(residue: int, sim: str, rep: int, chain: str) -> None:
    # Skip if analysis was already performed.
    if os.path.exists(f"ioncontacts/{residue}_{sim}_{rep}_{chain}.txt"):
        return

    base = f"../../sims/{sim}/{rep:02d}"
    u = MDAnalysis.Universe(f"{base}/CA.pdb", f"{base}/MD_conv.xtc")

    target = f"(chainID {chain} and resid {residue} and name {carboxylAtoms})"

    frames = []
    dists  = []

    num = getLambdaFileIndices(universe=u, resid=residue)[chains.index(chain)]
    data = loadxvg(f"{base}/cphmd-coord-{num}.xvg", dt=1000, b=0)

    lambdaTrajt = [val / 1000.0 for val in data[0]]  # debug
    lambdaTrajx = data[1]

    frame = 0
    for _ in u.trajectory:

        dist = distance_array(u.select_atoms(target), u.select_atoms("name NA")).min()
        frames.append(frame)
        dists.append(dist)
        frame += 1

    with open(f"ioncontacts/{residue}_{sim}_{rep}_{chain}.txt", 'w') as file:
        for idx in range(0, len(frames)):
            file.write(f"{frames[idx]} {dists[idx]} {lambdaTrajx[idx]}\n")

# Gather items:
items = []
for residue in residues:
    for sim in sims:
        for rep in reps:
            for chain in chains:
                items.append((residue, sim, rep, chain))

# Run multithreaded
pool = mp.Pool(processes=mp.cpu_count())
pool.starmap(task, items, chunksize=1)

# PLOTTING #####################################################################

for sim in sims:
    resTimeList = []

    for residue in residues:
        for rep in reps:
            for chain in chains:

                time = 0  # Occupancy time in ns.
                for mindist in loadCol(f"ioncontacts/{residue}_{sim}_{rep}_{chain}.txt", col=2):

                    if mindist < cutoff:
                        time += 1
                        continue
                    
                    if time != 0:
                        resTimeList.append(time)
                    
                    time = 0
            
            if time != 0:
                resTimeList.append(time)

    print(resTimeList)
    if resTimeList:
        print(f"{sim} {np.mean(resTimeList):.2f}±{stderr(resTimeList):.2f}")
