#!/usr/bin/env python3

import pickle
import MDAnalysis
from MDAnalysis.analysis.distances import distance_array
import multiprocessing as mp
import os

from science.utility import makeSuperDict

sims      = ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']
reps      = [1, 2, 3, 4]
chains    = ['A', 'B', 'C', 'D', 'E']

carboxylAtoms = 'name OE1 OE2 OD1 OD2 NE2 ND1'
polarAtoms    = 'name NZ NE2 OD2 OE1 OD1 NE OG1 OH ND2 OE2 OT1 NH2 O NH1 NE1 N ND1 OT2 OG'

# Speedup for dev purposes.
if not os.path.isfile('superIndices.obj'):

    # MULTIPROCESSING TASK.
    def task(sim, rep):
        contactFrameIndices = []

        p1 = f'../sims/{sim}/{rep:02d}/CA.pdb'
        p2 = f'../sims/{sim}/{rep:02d}/MD_conv.xtc'
        u  = MDAnalysis.Universe(p1, p2)

        for chain in chains:

            sel1 = u.select_atoms(f'chainID {chain} and resid 243 and {carboxylAtoms}')
            sel2 = u.select_atoms(f'chainID {chain} and resid 248 and {polarAtoms}')

            for ts in u.trajectory:
                if distance_array(sel1, sel2).min() < 4.0:
                    contactFrameIndices.append(ts.frame)

            # Dumping and loading using pickle is a hack to avoid using shared memory.
            pickle.dump(contactFrameIndices, open(f'{sim}_{rep}_{chain}.tmp', 'wb'))

    # PREPARE ITERABLES.
    items = []
    for sim in sims:
        for rep in reps:
            items.append((sim, rep))

    # RUN MULTITHREADED.
    pool = mp.Pool(processes=mp.cpu_count())
    pool.starmap(task, items, chunksize=1)

    # RECONSTRUCT.
    superIndices = makeSuperDict([sims, reps, chains, []])
    for sim in sims:
        for rep in reps:
            for chain in chains:
                superIndices[sim][rep][chain] = pickle.load(open(f'{sim}_{rep}_{chain}.tmp', 'rb'))
                os.remove(f'{sim}_{rep}_{chain}.tmp')

    pickle.dump(superIndices, open('superIndices.obj', 'wb'))

else:
    superIndices = pickle.load(open('superIndices.obj', 'rb'))
