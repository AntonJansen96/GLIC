#!/usr/bin/env python3

from science.parsing import loadxvg
from science.cphmd import getLambdaFileIndices
from science.cphmd import movingDeprotonation

import os
import MDAnalysis
import matplotlib
import matplotlib.pyplot as plt
import multiprocessing as mp

matplotlib.rcParams.update({'font.size': 20})

resids = [13, 14, 26, 31, 32, 35, 49, 55, 67, 69, 75, 82, 86, 88, 91, 97, 104, 115, 122, 127, 136, 145, 147, 153, 154, 161, 163, 177, 178, 181, 185, 222, 235, 243, 272, 277, 282]
# resids = [13]
sims   = ['6ZGD_7', '6ZGD_4', '4HFI_7', '4HFI_4']
fancy  = ["Closed, pH 7.0", "Closed, pH 4.0", "Open, pH 7.0", "Open, pH 4.0"]
reps   = [1, 2, 3, 4]
chains = ['A', 'B', 'C', 'D', 'E']

# DEFINE TASK
def task(resid: int, dummy: any):
    # if os.path.isfile(f"proto_{resid}.png"):
    #     return

    u = MDAnalysis.Universe("../../sims/4HFI_4/01/CA.pdb")

    nrows = 4
    ncols = 4

    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(17, 9), dpi=150)

    for row in range(0, nrows):
        for col in range(0, ncols):

            subplt = axs[row, col]

            store = []
            for idx in range(0, len(chains)):
                num = getLambdaFileIndices(universe=u, resid=resid)[idx]

                print(sims[col], f"{reps[row]:02d}", chains[idx], num)

                data = loadxvg(f"../../sims/{sims[col]}/{reps[row]:02d}/cphmd-coord-{num}.xvg", dt=1000, b=0)
                t = [val / 1000.0 for val in data[0]]
                
                if resid in [127, 235, 277]:  # histidine is reversed.
                    x = movingDeprotonation(data[1])
                else:
                    x = [1 - val for val in movingDeprotonation(data[1])]

                subplt.plot(t, x, linewidth=2)
                store.append(x)

            mins = []
            maxs = []
            for idx in range(0, len(t)):
                mins.append(min([store[0][idx], store[1][idx], store[2][idx], store[3][idx], store[4][idx]]))
                maxs.append(max([store[0][idx], store[1][idx], store[2][idx], store[3][idx], store[4][idx]]))
            subplt.fill_between(t, mins, maxs, color='lightblue', alpha=0.5)

            # Set x-lim and y-lim.
            subplt.set_ylim([-0.1, 1.1])
            subplt.set_xlim([0, 1000])

            # subplt.text(15, 4.0, str(row+1))

            # If we're not in the last row, do not show the xticks.
            if row != nrows - 1:
                subplt.set_xticks([])
            else:
                subplt.set_xlabel("Time (ns)")
                subplt.set_xticks([0, 250, 500, 750, 1000])

            # If we're not in the first column, do not show the yticks.
            if col != 0:
                subplt.set_yticks([])
            else:
                subplt.set_ylabel("Proto frac")
                subplt.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])

            # Add title to top row.
            if row == 0:
                subplt.set_title(fancy[col], size=20)

    fig.tight_layout(pad=0.2)
    fig.savefig(f"proto_{resid}.png")
    fig.clear()

# GATHER ITERABLES
items = []
for resid in resids:
    items.append((resid, 1.0))

# RUN MULTITHREADED
pool = mp.Pool(processes=mp.cpu_count())
pool.starmap(task, items, chunksize=1)
