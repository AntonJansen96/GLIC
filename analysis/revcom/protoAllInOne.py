#!/usr/bin/env python3

import os
import numpy as np
import MDAnalysis
import matplotlib
import matplotlib.pyplot as plt

from science.cphmd import getLambdaFileIndices, movingDeprotonation
from science.parsing import loadxvg, pickleDump, pickleLoad
from science.utility import makeSuperDict

matplotlib.rcParams.update({'font.size': 20})

def getFullResName(num: int) -> str:
    data = ['D13', 'E14', 'E26', 'D31', 'D32', 'E35', 'D49', 'D55', 'E67', 'E69', 'E75', 'E82', 'D86', 'D88', 'D91', 'D97', 'E104', 'D115', 'D122', 'H127', 'D136', 'D145', 'E147', 'D153', 'D154', 'D161', 'E163', 'E177', 'D178', 'E181', 'D185', 'E222', 'H235', 'E243', 'E272', 'H277', 'E282']
    for idx in range(0, len(data)):
        if num == int(data[idx][1:]):
            return data[idx]

def stderr(array):
    return np.std(array) / np.sqrt(len(array))

resids = [13, 14, 26, 31, 32, 35, 49, 55, 67, 69, 75, 82, 86, 88, 91, 97, 104, 115, 122, 127, 136, 145, 147, 153, 154, 161, 163, 177, 178, 181, 185, 222, 235, 243, 272, 277, 282, 282, 282, 282]
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

                    if resid in [127, 235, 277]:                                    # HIS
                        array = [1 - val for val in movingDeprotonation(data[1])]
                    else:                                                           # Rest
                        array = movingDeprotonation(data[1])

                    superData[resid][sim][rep][chain] = array

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
        
        comb = length * [0]
        sdev = length * [0]

        for idx in range(length):
            
            data = []
            for rep in reps:
                for chain in chains:
                    data.append(1 - superData[resid][sim][rep][chain][idx])

            comb[idx] = np.mean(data)
            sdev[idx] = stderr(data)

        ultraMean[resid][sim] = comb
        ultraSdev[resid][sim] = sdev

nrows = 10
ncols = 4

fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(17, 23), dpi=200)

row = 0
col = 0

for resid in resids:

    subplt = axs[row, col]

    # Plot the mean of twenty chains
    for sim in sims:
        x = ultraMean[resid][sim]
        t = range(len(x))

        lower = [0] * len(t)
        upper = [0] * len(t)

        for idx in t:
            lower[idx] = x[idx] - ultraSdev[resid][sim][idx]
            upper[idx] = x[idx] + ultraSdev[resid][sim][idx]

        if sim == '6ZGD_7':
            subplt.plot(t, x, linewidth=2, color='#8856a7')
            subplt.fill_between(t, lower, upper, color='#8856a7', alpha=0.5)

        if sim == '6ZGD_4':
            subplt.plot(t, x, linewidth=2, color='#9ebcda')
            subplt.fill_between(t, lower, upper, color='#9ebcda', alpha=0.5)

        if sim == '4HFI_7':
            subplt.plot(t, x, linewidth=2, color='#8856a7', linestyle='--')
            subplt.fill_between(t, lower, upper, color='#8856a7', alpha=0.5, hatch='//', edgecolor='w', lw=1)
    
        if sim == '4HFI_4':
            subplt.plot(t, x, linewidth=2, color='#9ebcda', linestyle='--')
            subplt.fill_between(t, lower, upper, color='#9ebcda', alpha=0.5, hatch='\\\\', edgecolor='w', lw=1)

    # Set x-lim and y-lim.
    subplt.set_ylim([-0.05, 1.05])
    subplt.set_xlim([0, 1000])

    # If we're not in the last row, do not show the xticks.
    subplt.set_xticks([200, 400, 600, 800])
    if row != nrows - 1:
        subplt.set_xticklabels([])
        subplt.xaxis.set_ticks_position('none')
    else:
        subplt.set_xlabel("Time (ns)")

    # If we're not in the first column, do not show the yticks.
    subplt.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    if col != 0:
        subplt.set_yticklabels([])
        subplt.yaxis.set_ticks_position('none')
    else:
        subplt.set_ylabel("Proto frac")

    # Add (horizontal) grid to all subplots.
    subplt.grid(True, linestyle='--', axis='y')

    # Add title to each subplot.
    subplt.set_title(getFullResName(resid), size=20)

    # Increment row and column indices.
    col += 1
    if col == ncols:
        row += 1
        col = 0

fig.tight_layout(pad=0.2)
fig.savefig(f"allinone.png")
os.system(f"convert allinone.png -trim allinone.png")
fig.clear()
