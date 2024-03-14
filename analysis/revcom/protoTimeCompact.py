#!/usr/bin/env python3

import os
import sys
import numpy as np
import MDAnalysis
import matplotlib
import matplotlib.pyplot as plt
import multiprocessing as mp

from science.cphmd import getLambdaFileIndices, movingDeprotonation
from science.parsing import loadxvg, pickleDump, pickleLoad
from science.utility import makeSuperDict

matplotlib.rcParams.update({'font.size': 20})  # Set global matplotlib font size.

# PARAMETERS

resids = [13, 14, 26, 31, 32, 35, 49, 55, 67, 69, 75, 82, 86, 88, 91, 97, 104, 115, 122, 127, 136, 145, 147, 153, 154, 161, 163, 177, 178, 181, 185, 222, 235, 243, 272, 277, 282]
sims   = ['6ZGD_7', '6ZGD_4', '4HFI_7', '4HFI_4']
fancy  = ["Closed, pH 7.0", "Closed, pH 4.0", "Open, pH 7.0", "Open, pH 4.0"]
reps   = [1, 2, 3, 4]
chains = ['A', 'B', 'C', 'D', 'E']

# GATHER ALL DATA INTO ONER SUPERDICT

if os.path.isfile("protoTimeCompact.obj"):
    
    superData = pickleLoad("protoTimeCompact.obj")
    print("Loaded pickled protoTimeCompact object...")

else:

    u = MDAnalysis.Universe("../../sims/4HFI_4/01/CA.pdb")

    # Holds moving deprotonations for all resids, sims, reps, chains.
    superData = makeSuperDict([resids, sims, reps, chains, []])  

    for resid in resids:
        for sim in sims:
            print(resid, sim)
            for rep in reps:
                for chain in chains:
                    num = getLambdaFileIndices(u, resid)[chains.index(chain)]
                    data = loadxvg(f"../../sims/{sim}/{rep:02d}/cphmd-coord-{num}.xvg", dt=1000)
                    superData[resid][sim][rep][chain] = movingDeprotonation(data[1])

    pickleDump(superData, "protoTimeCompact.obj")

# GET THE AVERAGE OVER THE RUNNING AVERAGES OF THE FIVE CHAINS

superMean = makeSuperDict([resids, sims, reps, []])

for resid in resids:
    for sim in sims:
        for rep in reps:

            length = len(superData[resid][sim][rep]['A'])
            comb = [0] * length

            for idx in range(length):
                a = 1 - superData[resid][sim][rep]['A'][idx]  # 1 - x because we go 
                b = 1 - superData[resid][sim][rep]['B'][idx]  # from deprotonation
                c = 1 - superData[resid][sim][rep]['C'][idx]  # to protonation!
                d = 1 - superData[resid][sim][rep]['D'][idx]
                e = 1 - superData[resid][sim][rep]['E'][idx]

                comb[idx] = np.mean([a, b, c, d, e])

            superMean[resid][sim][rep] = comb

# GET THE AVERAGE OVER THE AVERAGE OF THE FOUR REPLICATES

ultraMean = makeSuperDict([resids, sims, []])
ultraSdev = makeSuperDict([resids, sims, []])

for resid in resids:
    for sim in sims:
        
        l1 = len(superMean[resid][sim][1])
        l2 = len(superMean[resid][sim][2])
        l3 = len(superMean[resid][sim][3])
        l4 = len(superMean[resid][sim][4])
        length = min([l1, l2, l3, l4])
        comb = [0] * length
        sdev = [0] * length

        for idx in range(length):
            r1 = superMean[resid][sim][1][idx]
            r2 = superMean[resid][sim][2][idx]
            r3 = superMean[resid][sim][3][idx]
            r4 = superMean[resid][sim][4][idx]

            comb[idx] = np.mean([r1, r2, r3, r4])
            sdev[idx] = np.std([r1, r2, r3, r4])

        ultraMean[resid][sim] = comb
        ultraSdev[resid][sim] = sdev


resid = 13
sim = '6ZGD_7'

for rep in reps:
    t = range(len(superMean[resid][sim][rep]))
    x = superMean[resid][sim][rep]
    plt.plot(t, x, label=rep)

# Plot mean
t = range(len(ultraMean[resid][sim]))
x = ultraMean[resid][sim]
plt.plot(t, x, label='All', color='black')

# Plot sdev
lower = [0] * len(t)
upper = [0] * len(t)
for idx in t:
    lower[idx] = x[idx] - ultraSdev[resid][sim][idx]
    upper[idx] = x[idx] + ultraSdev[resid][sim][idx]
plt.fill_between(t, lower, upper, color='gray', alpha=0.3)

plt.legend()
plt.axis([0, 1000, -0.05, 1.05])
plt.savefig('test.png')    

sys.exit(0)







blocks = [[13, 14, 31, 32, 49, 55, 67, 69], [75, 86, 88, 91, 97, 104, 115, 122], [127, 136, 145, 147, 153, 154, 161, 163], [178, 181, 185, 222, 235, 272, 277, 282]]

for block in blocks:

    fig, axs = plt.subplots(nresids, nsims, figsize=(17, 9), dpi=200)








# DEFINE TASK
def task(resid: int, dummy: any):

    u = MDAnalysis.Universe("../../sims/4HFI_4/01/CA.pdb")

    nresids = len(reps) + 1  # 5, add a plus one for bottom row of 'All' panels.
    nsims = len(sims)          # 4

    fig, axs = plt.subplots(nresids, nsims, figsize=(17, 9), dpi=200)

    for simidx in range(0, nsims):

        # Stores the data for making the 'All' subplots.
        allValues = makeSuperDict([list(range(0, 1100)), []])

        for repidx in range(0, nresids):

            subplt = axs[repidx, simidx]

            if repidx < 4:

                store = []
                for idx in range(0, len(chains)):
                    num = getLambdaFileIndices(universe=u, resid=resid)[idx]

                    # Debug and user update.
                    print(sims[simidx], f"{reps[repidx]:02d}", chains[idx], num)

                    # Load data from associated cphmd-coord file.
                    data = loadxvg(f"../../sims/{sims[simidx]}/{reps[repidx]:02d}/cphmd-coord-{num}.xvg", dt=1000)
                    t = [val / 1000.0 for val in data[0]]

                    # Histidine is reversed.
                    if resid in [127, 235, 277]:
                        x = movingDeprotonation(data[1])
                    else:
                        x = [1 - val for val in movingDeprotonation(data[1])]

                    subplt.plot(t, x, linewidth=1.2, linestyle='--')
                    store.append(x)

                means = [0] * len(t)
                lower = [0] * len(t)
                upper = [0] * len(t)
                for idx in range(len(t)):
                    temp = [store[0][idx], store[1][idx], store[2][idx], store[3][idx], store[4][idx]]
                    mean = np.mean(temp)
                    sdev = np.std(temp)
                    means[idx] = mean
                    upper[idx] = mean + sdev
                    lower[idx] = mean - sdev

                    # Store data of the five chains for the 'All' plot.
                    allValues[idx] += temp  

                subplt.plot(t, means, color='black', linewidth=2.0)
                subplt.fill_between(t, lower, upper, color='gray', alpha=0.3)

            else:
                means = []
                lower = []
                upper = []
                finallength = 0

                for idx in range(0, 1100):
                    values = allValues[idx]

                    if len(values) < 10:    # This is because not all sims ran equally long.
                        finallength = idx   # Some stopped slightly before 1000 ns.
                        break

                    mean = np.mean(values)
                    sdev = np.std(values)
                    means.append(mean)
                    lower.append(mean - sdev)
                    upper.append(mean + sdev)

                subplt.plot(range(finallength), means, color='black', linewidth=2.0)
                subplt.fill_between(range(finallength), lower, upper, color='gray', alpha=0.3)

            # Set x-lim and y-lim.
            subplt.set_ylim([-0.05, 1.05])
            subplt.set_xlim([0, 1000])

            if repidx < 4:
                subplt.text(15, 0.85, str(repidx+1))
            else:
                subplt.text(15, 0.85, 'All')

            # If we're not in the last row, do not show the xticks.
            subplt.set_xticks([200, 400, 600, 800])
            if repidx != nresids - 1:
                subplt.set_xticklabels([])
                subplt.xaxis.set_ticks_position('none')
            else:
                subplt.set_xlabel("Time (ns)")

            # If we're not in the first column, do not show the yticks.
            subplt.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
            if simidx != 0:
                subplt.set_yticklabels([])
                subplt.yaxis.set_ticks_position('none')
            else:
                subplt.set_ylabel("Proto frac")

            # Add (horizontal) grid to all subplots.
            subplt.grid(True, linestyle='--', axis='y')

            # Add title to top row.
            if repidx == 0:
                subplt.set_title(fancy[simidx], size=20)

    fig.tight_layout(pad=0.2)
    fig.savefig(f"proto_{resid}.png")
    os.system(f"convert proto_{resid}.png -trim proto_{resid}.png")
    fig.clear()

# GATHER ITERABLES
items = []
for resid in resids:
    items.append((resid, 1.0))

# RUN MULTITHREADED
pool = mp.Pool(processes=mp.cpu_count())
pool.starmap(task, items, chunksize=1)
