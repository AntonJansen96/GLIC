#!/bin/python3

import matplotlib
import matplotlib.pyplot as plt

from science.parsing import loadxvg
from science.parsing import loadCol

# Set global font size for figures.
matplotlib.rcParams.update({'font.size': 14})

cutoff = 0.40  # nm

for target in [35]:

    for rep in [1]:

        #? IN THIS LOOP WE MAKE ONE PANEL
        for sim in ['4HFI_4']:
        # for sim in ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']:
            # for chain in ['A']:
            for chain in ['A', 'B', 'C', 'D', 'E']:

                #? IN THIS LOOP WE MAKE ONE PLOT

                # Load the top 6 contacts for said residue from the contacts/target.dat file
                topList = loadCol('contacts/{}.dat'.format(target), col=1, header=1)[:6]
                # If there is a Na+, remove it. If not, remove the last element.
                if 'Na+' in topList:
                    topList.remove('Na+')
                else:
                    topList = topList[:5]
                topList.append('Na+')

                A = []  # This is for y-ticks
                B = []  # This is for y-ticks
                for ii in range(0, len(topList)):

                    if topList[ii] == 'Na+':
                        continue  # continue here as other program not currently finished yet
                        fname = 'time/{}_{}_{}_{}_NA.xvg'.format(sim, rep, target, chain)
                        print(fname)

                    elif topList[ii][-1] in ['c', 'p']:
                        fname = 'time/{}_{}_{}_{}_{}.xvg'.format(sim, rep, target, chain, int(topList[ii][1:-1]))
                        print(fname)

                    else:
                        fname = 'time/{}_{}_{}_{}_{}.xvg'.format(sim, rep, target, chain, int(topList[ii][1:]))
                        print(fname)

                    D = loadxvg(fname=fname)
                    t = [val / 1000.0 for val in D[0]]  # ps -> ns
                    x = D[1]

                    for jj in range(0, len(x)):
                        if x[jj] < cutoff:
                            x[jj] = ii + 1
                        else:
                            x[jj] = -1

                    plt.scatter(t, x, marker='_')
                    A.append(ii + 1)
                    B.append(topList[ii])

                plt.axis([0, 1000, 0, 6])
                plt.xlabel('Time (ns)')
                plt.yticks(A, B)
                plt.title('{} in chain {}'.format(target, chain))
                plt.tight_layout()
                plt.savefig('time/{}_{}.png'.format(sim, chain))
                plt.clf()
                plt.close()
