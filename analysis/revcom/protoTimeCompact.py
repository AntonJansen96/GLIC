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
    """Returns 'D13' if num=13, etc."""
    data = ['D13', 'E14', 'E26', 'D31', 'D32', 'E35', 'D49', 'D55', 'E67', 'E69', 'E75', 'E82', 'D86', 'D88', 'D91', 'D97', 'E104', 'D115', 'D122', 'H127', 'D136', 'D145', 'E147', 'D153', 'D154', 'D161', 'E163', 'E177', 'D178', 'E181', 'D185', 'E222', 'H235', 'E243', 'E272', 'H277', 'E282']
    for idx in range(0, len(data)):
        if num == int(data[idx][1:]):
            return data[idx]

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

blocks = [[13, 14, 31, 32, 49, 55, 67, 69], [75, 86, 88, 91, 97, 104, 115, 122], [127, 136, 145, 147, 153, 154, 161, 163], [178, 181, 185, 222, 235, 272, 277, 282]]

for block in blocks:

    nresids = len(block)
    nsims = 4

    fig, axs = plt.subplots(nresids, nsims, figsize=(17, 18), dpi=200)

    for simidx in range(0, nsims):
        
        for residx in range(0, nresids):

            subplt = axs[residx, simidx]

            # Plot the four replicas
            for rep in reps:
                x = superMean[block[residx]][sims[simidx]][rep]
                t = range(len(x))
                subplt.plot(t, x, linewidth=1.2, linestyle='--')

            # Plot the mean of the four replicas.
            x = ultraMean[block[residx]][sims[simidx]]
            t = range(len(x))
            subplt.plot(t, x, linewidth=1.5, color='black')

            # Plot the standard deviation of the four replicas.
            lower = [0] * len(t)
            upper = [0] * len(t)
            for idx in t:
                lower[idx] = x[idx] - ultraSdev[block[residx]][sims[simidx]][idx]
                upper[idx] = x[idx] + ultraSdev[block[residx]][sims[simidx]][idx]
            subplt.fill_between(t, lower, upper, color='gray', alpha=0.3)

            # Set x-lim and y-lim.
            subplt.set_ylim([-0.05, 1.05])
            subplt.set_xlim([0, 1000])

            # If we're not in the last row, do not show the xticks.
            subplt.set_xticks([200, 400, 600, 800])
            if residx != nresids - 1:
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
                subplt.set_ylabel(f"{getFullResName(block[residx])}", size=24)
                # subplt.text(-400, 0.44)

            # Add (horizontal) grid to all subplots.
            subplt.grid(True, linestyle='--', axis='y')

            # Add title to top row.
            if residx == 0:
                subplt.set_title(fancy[simidx], size=20)

    fig.tight_layout(pad=0.2)
    fig.savefig(f"compact_{block[0]}.png")
    os.system(f"convert compact_{block[0]}.png -trim compact_{block[0]}.png")
    fig.clear()
