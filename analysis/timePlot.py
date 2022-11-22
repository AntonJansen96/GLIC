#!/bin/python3

import MDAnalysis
import matplotlib.pyplot as plt

from science.parsing import loadxvg
from science.parsing import loadCol
from science.cphmd import getLambdaFileIndices
from science.utility import triplet2letter

# Set global font size for figures.
# import matplotlib
# matplotlib.rcParams.update({'font.size': 14})

cutoff = 0.40  # contact distance cutoff (nm)
protoC = 0.2   # protonation binning cutoff


def calcLetter(letter, type):
    N = ['A', 'B', 'C', 'D', 'E']
    C = ['E', 'A', 'B', 'C', 'D']
    P = ['B', 'C', 'D', 'E', 'A']

    if type == 'c':
        for idx in range(0, len(N)):
            if letter == N[idx]:
                return C[idx]

    elif type == 'p':
        for idx in range(0, len(N)):
            if letter == N[idx]:
                return P[idx]


for target in [13, 14, 26, 31, 32, 35, 49, 55, 67, 69, 75, 82, 86, 88, 91, 97, 104, 115, 122, 127, 136, 145, 147, 153, 154, 161, 163, 177, 178, 181, 185, 222, 235, 243, 272, 277, 282]:

    # Put this in outer loop for speedup
    lambdaIndices = getLambdaFileIndices('../sims/4HFI_4/01/CA.pdb', target)

    for rep in [1, 2, 3, 4]:

        #? IN THIS LOOP WE MAKE ONE SUPERPLOT

        fig, axs = plt.subplots(nrows=5, ncols=4, figsize=(17, 11), dpi=400)

        sims = ['6ZGD_7', '6ZGD_4', '4HFI_7', '4HFI_4']
        for sim in sims:

            chain1 = ['A', 'B', 'C', 'D', 'E']
            for chain in chain1:

                #? IN THIS LOOP WE MAKE ONE SUBPLOT

                subplt = axs[chain1.index(chain), sims.index(sim)]

                A = []  # This is for y-ticks (positions)
                B = []  # This is for y-ticks (labels)

                #! Do the protonation part (bottom-most line).

                lambdaIndex = lambdaIndices[chain1.index(chain)]
                print('../sims/{}/{:02d}/cphmd-coord-{}.xvg'.format(sim, rep, lambdaIndex))  # debug

                D = loadxvg('../sims/{}/{:02d}/cphmd-coord-{}.xvg'.format(sim, rep, lambdaIndex), dt=1000, b=0)
                t = [val / 1000.0 for val in D[0]]  # ps -> ns
                x = D[1]

                # Do the binning
                for jj in range(0, len(x)):
                    if x[jj] < protoC:
                        x[jj] = 1
                    else:
                        x[jj] = -1

                subplt.scatter(t, x, marker='|', linewidth=1.0)
                A.append(1)
                B.append('proto')

                #! Do the residues/Na+ part.

                # Load the top 6 contacts for said residue from the contacts/target.dat file
                topList = loadCol('contacts/{}.dat'.format(target), col=1, header=1)[:6]
                # If there is a Na+, remove it. If not, remove the last element.
                # Then, include Na+ by default at the end.
                if 'Na+' in topList:
                    topList.remove('Na+')
                else:
                    topList = topList[:5]
                topList.append('Na+')

                for ii in range(0, len(topList)):

                    if topList[ii] == 'Na+':  # SODIUM
                        fname = 'time/{}_{}_{}{}_NA.xvg'.format(sim, rep, target, chain)

                    elif topList[ii][-1] == 'c':  # INTER-CHAIN CASE C
                        xx = calcLetter(chain, type='c')
                        fname = 'time/{}_{}_{}{}_{}{}.xvg'.format(sim, rep, target, chain, int(topList[ii][1:-1]), xx)

                    elif topList[ii][-1] == 'p':  # INTER-CHAIN CASE P
                        xx = calcLetter(chain, type='p')
                        fname = 'time/{}_{}_{}{}_{}{}.xvg'.format(sim, rep, target, chain, int(topList[ii][1:-1]), xx)

                    else:  # SAME CHAIN
                        fname = 'time/{}_{}_{}{}_{}{}.xvg'.format(sim, rep, target, chain, int(topList[ii][1:]), chain)

                    print(fname)  # debug
                    D = loadxvg(fname=fname)
                    t = [val / 1000.0 for val in D[0]]  # ps -> ns
                    x = D[1]

                    # Do the binning
                    for jj in range(0, len(x)):
                        if x[jj] < cutoff:
                            x[jj] = ii + 2     # The first one should be one
                        else:                  # higher than protonation.
                            x[jj] = -1

                    subplt.scatter(t, x, marker='|', linewidth=1.0)
                    A.append(ii + 2)
                    B.append(topList[ii])

                subplt.axis([0, 1000, 0, 8])
                subplt.set_yticks(A, B)

        #! Finish the panel

        # Set outer xticks (top) (closed, pH 7, etc.)
        simNames = ['closed, pH 7', 'closed, pH 4', 'open, pH 7', 'open, pH 4']
        for jj in range(0, len(sims)):
            axs[0, jj].set_title(simNames[jj], size=16, pad=10)

        # Set outer xticks (bottom) (Time (ns))
        for jj in range(0, len(sims)):
            axs[4, jj].set_xlabel('Time (ns)', size=14)

        # Set outer yticks (A, B, C, D, E)
        for ii in range(0, len(chain1)):
            axs[ii, 0].set_ylabel(chain1[ii], rotation=0, size=16, labelpad=20)

        # Set the super title (get the restype letter as well)
        u = MDAnalysis.Universe('../sims/4HFI_4/01/CA.pdb')
        triplet = u.select_atoms('chainID A and resid {}'.format(target)).residues.resnames[0]
        letter = triplet2letter(triplet)
        fig.suptitle('Residue {}{} (replica {})'.format(letter, target, rep), size=20, y=0.995)

        # Save and clear
        fig.tight_layout(pad=0.6)
        fig.savefig('time/{}_{}.png'.format(target, rep))
        fig.clf()
