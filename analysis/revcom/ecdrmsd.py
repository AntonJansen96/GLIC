#!/usr/bin/env python3

import os
import MDAnalysis
import MDAnalysis.analysis.rms
import multiprocessing as mp
import matplotlib
import matplotlib.pyplot as plt

from science.parsing import loadCol

# sims = ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']
sims = ['6ZGD_7', '6ZGD_4', '4HFI_7', '4HFI_4']
fancy = ["Closed, pH 7.0", "Closed, pH 4.0", "Open, pH 7.0", "Open, pH 4.0"]
reps = [1, 2, 3, 4]
chains = ['A', 'B', 'C', 'D', 'E']

# sel = "resid 18 to 200 and name CA"
sel = "((resid 15 to 48) or (resid 66 to 192)) and name CA"  # Urska21, checked in VMD.

# DEFINE TASK
def task(sim, rep, chain):
    if not os.path.exists(f"{sim}_{rep}_{chain}.txt"):
        path1 = f"../../sims/{sim}/{rep:02d}/CA.pdb"
        path2 = f"../../sims/{sim}/{rep:02d}/MD_conv.xtc"
        u = MDAnalysis.Universe(path1, path2)
        R = MDAnalysis.analysis.rms.RMSD(u, select=f"chainID {chain} and ({sel})")
        R.run()

        t = [val / 1000.0 for val in R.rmsd.T[1]]
        x = R.rmsd.T[2]
        with open(f"{sim}_{rep}_{chain}.txt", 'w+') as file:
            for idx in range(0, len(t)):
                file.write(f"{t[idx]} {x[idx]}\n")

# GATHER ITERABLES
items = []
for sim in sims:
    for rep in reps:
        for chain in chains:
            items.append((sim, rep, chain))

# RUN MULTITHREADED
pool = mp.Pool(processes=mp.cpu_count())
pool.starmap(task, items, chunksize=1)

# DO PLOTTING
matplotlib.rcParams.update({'font.size': 22})

nrows = 4
ncols = 4

fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(17, 11), dpi=200)

for row in range(0, nrows):
    for col in range(0, ncols):

        subplt = axs[row, col]

        for chain in chains:
            t = loadCol(f"{sims[col]}_{reps[row]}_{chain}.txt", col=1)
            x = loadCol(f"{sims[col]}_{reps[row]}_{chain}.txt", col=2)
            subplt.plot(t, x, linewidth=0.5)

        # Set x-lim and y-lim.
        subplt.set_ylim([0.5, 4.5])
        subplt.set_xlim([0, 1000])

        subplt.text(15, 4.0, str(row+1))

        # If we're not in the last row, do not show the xticks.
        if row != nrows - 1:
            subplt.set_xticks([])
        else:
            subplt.set_xlabel("Time (ns)")
            subplt.set_xticks([200, 400, 600, 800])

        # If we're not in the first column, do not show the yticks.
        if col != 0:
            subplt.set_yticks([])
        else:
            subplt.set_ylabel(r'ECD RMSD ($\AA$)', size=20)
            subplt.set_yticks([1, 2, 3, 4])

        # Add title to top row.
        if row == 0:
            subplt.set_title(fancy[col])

fig.tight_layout(pad=0.3)
fig.savefig('ecdrmsd.png')
os.system("convert ecdrmsd.png -trim ecdrmsd.png")
