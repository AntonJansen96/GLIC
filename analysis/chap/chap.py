#!/bin/python3

import json
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# Set global font size for figures.
matplotlib.rcParams.update({'font.size': 14})


def doPlot1(sims, reps, title, fname):

    plt.figure(figsize=(5, 7))

    for sim in sims:
        for rep in reps:
            # Load protein.json data file.
            data = json.load(open('{}/{:02d}/protein.json'.format(sim, rep)))

            # Plot time-averaged pore radius.
            r = [i * 10 for i in data['pathwayProfile']['radiusMean']]
            z = [i * 10 for i in data['pathwayProfile']['s']]
            plt.plot(r, z, lw=2, label='{}_{}'.format(sim, rep))

            # Plot associated standard deviation as a shaded region.
            s = [i * 10 for i in data['pathwayProfile']['radiusSd']]
            A = np.array(z) + np.array(s)
            B = np.array(z) - np.array(s)
            plt.fill_between(r, A, B, facecolor='#000000', alpha=0.2)

    plt.xlabel('Radius (Å)')
    plt.ylabel('Z-coordinate (Å)')
    plt.axis([0, 6, -20, 20])
    plt.legend(fontsize='small')
    plt.title(title)
    plt.tight_layout()
    plt.savefig(fname)
    plt.clf()
    plt.close()


sims = ['6ZGD_7', '4HFI_7', '6ZGD_4', '4HFI_4']
reps = [1, 2, 3, 4]

for sim in sims:
    doPlot1(['{}'.format(sim)], reps, 'Simulation {}'.format(sim), 'profile_{}.png'.format(sim))

for rep in reps:
    doPlot1(sims, [rep], 'Replica {}'.format(rep), 'profile_{}.png'.format(rep))

#? DO THE STANDARD DEVATION/ERROR INSTEAD BY AVERAGING OVER FOUR REPLICAS: #####
#* Fuck this shit: z-lists are not aligned for different replicas.

# from science.utility import makeSuperDict

# sims = ['4HFI_4']
# reps = [1, 2, 3, 4]

# # Create and fill data structure for the trajectories.
# radiusData = makeSuperDict([sims, reps, []])
# for sim in sims:
#     for rep in reps:
#         data = json.load(open('{}/{:02d}/protein.json'.format(sim, rep)))
#         radiusData[sim][rep] = [i * 10 for i in data['pathwayProfile']['radiusMean']]

# plt.figure(figsize=(5, 7))

# # Loop through the simulations
# for sim in sims:

#     meanList = []
#     stdrList = []

#     # Loop through the values
#     for idx in range(0, len(radiusData[sim][1])):

#         rep1 = radiusData[sim][1][idx]
#         rep2 = radiusData[sim][2][idx]
#         rep3 = radiusData[sim][3][idx]
#         rep4 = radiusData[sim][4][idx]

#         meanList.append(np.mean([rep1, rep2, rep3, rep4]))
#         stdrList.append(np.std( [rep1, rep2, rep3, rep4]))

#     # Plot replica-average pore radius.
#     data = json.load(open('{}/{:02d}/protein.json'.format(sim, 1)))
#     z    = [i * 10 for i in data['pathwayProfile']['s']]
#     plt.plot(meanList, z, lw=2, label='{}'.format(sim))

#     # Plot associated standard deviation as a shaded region.
#     A = np.array(z) + np.array(stdrList)
#     B = np.array(z) - np.array(stdrList)
#     plt.fill_between(meanList, A, B, facecolor='#000000', alpha=0.2)

# plt.xlabel('Radius (Å)')
# plt.ylabel('Z-coordinate (Å)')
# plt.axis([0, 6, -20, 20])
# plt.legend(fontsize='small')
# plt.title('custom')
# plt.tight_layout()
# plt.savefig('profile2.png')
# plt.clf()
# plt.close()


# data1 = json.load(open('{}/{:02d}/protein.json'.format('4HFI_4', 1)))
# z1 = [i * 10 for i in data1['pathwayProfile']['s']]
# data2 = json.load(open('{}/{:02d}/protein.json'.format('4HFI_4', 2)))
# z2 = [i * 10 for i in data2['pathwayProfile']['s']]
# data3 = json.load(open('{}/{:02d}/protein.json'.format('4HFI_4', 3)))
# z3 = [i * 10 for i in data3['pathwayProfile']['s']]
# data4 = json.load(open('{}/{:02d}/protein.json'.format('4HFI_4', 4)))
# z4 = [i * 10 for i in data4['pathwayProfile']['s']]

# try:
#     for idx in range(0, 10000000 + 1):
#         print(z1[idx], z2[idx], z3[idx], z4[idx], z4[idx+1]-z4[idx])
# except IndexError:
#     pass
