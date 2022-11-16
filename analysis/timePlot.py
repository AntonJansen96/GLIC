#!/bin/python3

import matplotlib
import matplotlib.pyplot as plt

from science.parsing import loadxvg
from science.parsing import loadCol
from science.cphmd import getLambdaFileIndices
from science.cphmd import protonation

# Set global font size for figures.
matplotlib.rcParams.update({'font.size': 14})

cutoff = 0.40  # contact distance cutoff (nm)
protoC = 0.2   # protonation binning cutoff

for target in [35]:
    for rep in [1]:

        #? IN THIS LOOP WE MAKE ONE PANEL

        for sim in ['6ZGD_7']:

            lambdaIndices = getLambdaFileIndices('../sims/{}/{:02d}/CA.pdb'.format(sim, rep), target)
            chain1 = ['A', 'B']

            for chain in chain1:

                #? IN THIS LOOP WE MAKE ONE PLOT

                A = []  # This is for y-ticks (positions)
                B = []  # This is for y-ticks (labels)

                #! Do the protonation part (bottom-most line).

                lambdaIndex = lambdaIndices[chain1.index(chain)]
                print('../sims/{}/{:02d}/cphmd-coord-{}.xvg'.format(sim, rep, lambdaIndex))  # debug

                D = loadxvg('../sims/{}/{:02d}/cphmd-coord-{}.xvg'.format(sim, rep, lambdaIndex), dt=1000, b=0)
                t = [val / 1000.0 for val in D[0]]  # ps -> ns
                x = D[1]
                print(protonation(x))  # debug

                # Do the binning
                for jj in range(0, len(x)):
                    if x[jj] < protoC:
                        x[jj] = 1
                    else:
                        x[jj] = -1

                plt.scatter(t, x, marker='|', linewidth=1.0)
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

                    if topList[ii] == 'Na+':
                        fname = 'time/{}_{}_{}_{}_NA.xvg'.format(sim, rep, target, chain)
                        print(fname)  # debug

                    elif topList[ii][-1] in ['c', 'p']:
                        fname = 'time/{}_{}_{}_{}_{}.xvg'.format(sim, rep, target, chain, int(topList[ii][1:-1]))
                        print(fname)  # debug

                    else:
                        fname = 'time/{}_{}_{}_{}_{}.xvg'.format(sim, rep, target, chain, int(topList[ii][1:]))
                        print(fname)  # debug

                    D = loadxvg(fname=fname)
                    t = [val / 1000.0 for val in D[0]]  # ps -> ns
                    x = D[1]

                    # Do the binning
                    for jj in range(0, len(x)):
                        if x[jj] < cutoff:
                            x[jj] = ii + 2
                        else:
                            x[jj] = -1

                    plt.scatter(t, x, marker='|', linewidth=1.0)
                    A.append(ii + 2)
                    B.append(topList[ii])

                plt.axis([0, 1000, 0, 8])
                plt.xlabel('Time (ns)')
                plt.yticks(A, B)
                plt.title('{} in chain {}'.format(target, chain))
                plt.tight_layout()
                plt.savefig('time/{}_{}.png'.format(sim, chain))
                plt.clf()
                plt.close()
