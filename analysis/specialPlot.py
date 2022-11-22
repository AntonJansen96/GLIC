#!/bin/python3

import matplotlib
import matplotlib.pyplot as plt
import copy
import numpy as np
import MDAnalysis

from science.parsing import loadCol
from science.parsing import loadxvg
from science.cphmd import getLambdaFileIndices
from science.utility import triplet2letter

# Set global font size for figures.
matplotlib.rcParams.update({'font.size': 14})


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


def stderr(array):
    return np.std(array) / np.sqrt(len(array))


cutoff  = 0.40  # contact distance cutoff (nm)
protoC  = 0.2   # protonation binning cutoff

sims    = ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']
reps    = [1, 2, 3, 4]
chains  = ['A', 'B', 'C', 'D', 'E']
targets = [13, 14, 26, 31, 32, 35, 49, 55, 67, 69, 75, 82, 86, 88, 91, 97, 104, 115, 122, 127, 136, 145, 147, 153, 154, 161, 163, 177, 178, 181, 185, 222, 235, 243, 272, 277, 282]

for target in targets:

    #? GET LAMBDAINDICES AND TOPLIST
    #? lambdaIndices are the same for all. topList contains the contact names.

    lambdaIndices = getLambdaFileIndices('../sims/4HFI_4/01/CA.pdb', target)

    # Load the top 6 contacts for said residue from the contacts/target.dat file
    topList = loadCol('contacts/{}.dat'.format(target), col=1, header=1)[:6]
    # If there is a Na+, remove it. If not, remove the last element.
    # Then, include Na+ by default at the end.
    if 'Na+' in topList:
        topList.remove('Na+')
    else:
        topList = topList[:5]
    topList.append('Na+')

    #? CONSTRUCT THE SUPERDATA
    #? (dict of sims of reps of chains of topResidues of (mindist) lists)

    A = {}
    for key in topList:
        A[key] = []
    B = {}
    for key in chains:
        B[key] = copy.deepcopy(A)
    C = {}
    for key in reps:
        C[key] = copy.deepcopy(B)
    superData = {}
    for key in sims:
        superData[key] = copy.deepcopy(C)

    # print(superData)  # debug

    #? CONSTRUCT THE SUPERLAMBDA
    #? (dict of sims of reps of chains of lists of coordinates)

    A = {}
    for key in chains:
        A[key] = []
    B = {}
    for key in reps:
        B[key] = copy.deepcopy(A)
    superLambda = {}
    for key in sims:
        superLambda[key] = copy.deepcopy(B)

    # print(superLambda)  # debug

    #? GATHER ALL THE DATA INTO SUPERDATA AND SUPERLAMBDA

    for sim in sims:
        for rep in reps:
            for chain in chains:

                #! Fill superLambda
                fname = '../sims/{}/{:02d}/cphmd-coord-{}.xvg'.format(sim, rep, lambdaIndices[chains.index(chain)])
                superLambda[sim][rep][chain] = loadxvg(fname, dt=1000)[1]
                print(fname, "(chain = {})".format(chain))  # debug

                #! Fill superData
                for name in topList:
                    if name == 'Na+':  # SODIUM
                        fname = 'time/{}_{}_{}{}_NA.xvg'.format(sim, rep, target, chain)

                    elif name[-1] == 'c':  # INTER-CHAIN CASE C
                        xx = calcLetter(chain, type='c')
                        fname = 'time/{}_{}_{}{}_{}{}.xvg'.format(sim, rep, target, chain, int(name[1:-1]), xx)

                    elif name[-1] == 'p':  # INTER-CHAIN CASE P
                        xx = calcLetter(chain, type='p')
                        fname = 'time/{}_{}_{}{}_{}{}.xvg'.format(sim, rep, target, chain, int(name[1:-1]), xx)

                    else:  # SAME CHAIN
                        fname = 'time/{}_{}_{}{}_{}{}.xvg'.format(sim, rep, target, chain, int(name[1:]), chain)

                    superData[sim][rep][chain][name] = loadxvg(fname)[1]
                    print(fname)  # debug

    # print(superLambda)  # debug
    # print(superData)    # debug

    #! DO THE PH MIXING PART

    bins = ['deproto_closed', 'proto_closed', 'deproto_open', 'proto_open']

    #? CONSTRUCT THE SUPERRESULT
    # (dict of bins of topresidues of (empty) meanLists)
    A = {}
    for key in topList:
        A[key] = []
    superResult = {}
    for key in bins:
        superResult[key] = copy.deepcopy(A)

    # print(superResult)  # debug

    #? FILL THE SUPERRESULT

    for bin in bins:
        for name in topList:
            for rep in reps:
                for chain in chains:

                    frameCount = 0  # Number of frames where we have the right protonation state.
                    eventCount = 0  # Number of frames where we have, in addition, a contact.

                    # use 998 here instead of 1001 as not all simulations were
                    # run until exactly 1000ns.

                    if bin == 'deproto_closed':
                        for sim in ['6ZGD_4', '6ZGD_7']:  # If closed...
                            for idx in range(0, 998):
                                Lx = superLambda[sim][rep][chain][idx]
                                Dx   = superData[sim][rep][chain][name][idx]

                                if Lx > 1 - protoC:  # If deprotonated...
                                    frameCount += 1
                                    if Dx < cutoff:  # If within cutoff...
                                        eventCount += 1

                    elif bin == 'proto_closed':
                        for sim in ['6ZGD_4', '6ZGD_7']:  # If closed...
                            for idx in range(0, 998):
                                Lx = superLambda[sim][rep][chain][idx]
                                Dx   = superData[sim][rep][chain][name][idx]

                                if Lx < protoC:  # If protonated...
                                    frameCount += 1
                                    if Dx < cutoff:  # If within cutoff...
                                        eventCount += 1

                    elif bin == 'deproto_open':
                        for sim in ['4HFI_4', '4HFI_7']:  # If open...
                            for idx in range(0, 998):
                                Lx = superLambda[sim][rep][chain][idx]
                                Dx   = superData[sim][rep][chain][name][idx]

                                if Lx > 1 - protoC:  # If deprotonated...
                                    frameCount += 1
                                    if Dx < cutoff:  # If within cutoff...
                                        eventCount += 1

                    elif bin == 'proto_open':
                        for sim in ['4HFI_4', '4HFI_7']:  # If open...
                            for idx in range(0, 998):
                                Lx = superLambda[sim][rep][chain][idx]
                                Dx   = superData[sim][rep][chain][name][idx]

                                if Lx < protoC:  # If protonated...
                                    frameCount += 1
                                    if Dx < cutoff:  # If within cutoff...
                                        eventCount += 1

                    if frameCount == 0:
                        frac = 0
                    else:
                        frac = eventCount / float(frameCount)

                    superResult[bin][name].append(frac)

    #? MAKE THE BARPLOT

    # Initialize required data structures
    meanList = {}
    serrList = {}
    for key in bins:
        meanList[key] = []
        serrList[key] = []

    # Fill required data structures
    for bin in bins:
        for name in topList:
            meanList[bin].append(np.mean(superResult[bin][name]))
            serrList[bin].append(stderr(superResult[bin][name]))

    # print(meanList)  # debug

    # Some plotting stuff
    width = 0.2
    x     = np.arange(len(topList))
    fig   = plt.figure(figsize=(10, 6))
    ax    = fig.add_subplot()

    # deproto_closed
    mean4 = meanList['deproto_closed']
    serr4 = serrList['deproto_closed']
    ax.bar(     x - width * 1.5, mean4, width, color='#8856a7', label='closed, deprotonated')
    ax.errorbar(x - width * 1.5, mean4, serr4, color='#8856a7', fmt='none', capsize=6, linewidth=2)

    # proto_closed
    mean3 = meanList['proto_closed']
    serr3 = serrList['proto_closed']
    ax.bar(     x - width / 2.0, mean3, width, color='#9ebcda', label='closed, protonated')
    ax.errorbar(x - width / 2.0, mean3, serr3, color='#9ebcda', fmt='none', capsize=6, linewidth=2)

    # deproto_open
    mean2 = meanList['deproto_open']
    serr2 = serrList['deproto_open']
    ax.bar(     x + width / 2.0, mean2, width, color='#8856a7', label='open, deprotonated', edgecolor='w', lw=1, hatch='//')
    ax.errorbar(x + width / 2.0, mean2, serr2, color='#8856a7', fmt='none', capsize=6, linewidth=2)

    # proto_open
    mean1 = meanList['proto_open']
    serr1 = serrList['proto_open']
    ax.bar(     x + width * 1.5, mean1, width, color='#9ebcda', label='open, protonated', edgecolor='w', lw=1, hatch='\\\\')
    ax.errorbar(x + width * 1.5, mean1, serr1, color='#9ebcda', fmt='none', capsize=6, linewidth=2)

    ax.set_xticks(x, topList)
    ax.legend(loc=1, prop={'size': 12})

    # Do the title
    u = MDAnalysis.Universe('../sims/4HFI_4/01/CA.pdb')
    triplet = u.select_atoms('chainID A and resid {}'.format(target)).residues.resnames[0]
    letter = triplet2letter(triplet)
    plt.title('Residue {}{} (pH mixing)'.format(letter, target))

    plt.ylim(0, 1.1)
    plt.ylabel('Contact occupancy')
    plt.tight_layout()
    plt.savefig('contacts/{}_phmix.png'.format(target))
    plt.clf()
    plt.close()

    #! NOW WE DO THE STRUCTURE MIXING PART

    bins = ['deproto_ph7', 'proto_ph7', 'deproto_ph4', 'proto_ph4']

    #? CONSTRUCT THE SUPERRESULT
    # (dict of bins of topresidues of (empty) meanLists)
    A = {}
    for key in topList:
        A[key] = []
    superResult = {}
    for key in bins:
        superResult[key] = copy.deepcopy(A)

    # print(superResult)  # debug

    #? FILL THE SUPERRESULT

    for bin in bins:
        for name in topList:
            for rep in reps:
                for chain in chains:

                    frameCount = 0  # Number of frames where we have the right protonation state.
                    eventCount = 0  # Number of frames where we have, in addition, a contact.

                    # use 998 here instead of 1001 as not all simulations were
                    # run until exactly 1000ns.

                    if bin == 'deproto_ph7':
                        for sim in ['4HFI_7', '6ZGD_7']:  # If pH 7...
                            for idx in range(0, 998):
                                Lx = superLambda[sim][rep][chain][idx]
                                Dx   = superData[sim][rep][chain][name][idx]

                                if Lx > 1 - protoC:  # If deprotonated...
                                    frameCount += 1
                                    if Dx < cutoff:  # If within cutoff...
                                        eventCount += 1

                    elif bin == 'proto_ph7':
                        for sim in ['4HFI_7', '6ZGD_7']:  # If pH 7...
                            for idx in range(0, 998):
                                Lx = superLambda[sim][rep][chain][idx]
                                Dx   = superData[sim][rep][chain][name][idx]

                                if Lx < protoC:  # If protonated...
                                    frameCount += 1
                                    if Dx < cutoff:  # If within cutoff...
                                        eventCount += 1

                    elif bin == 'deproto_ph4':
                        for sim in ['4HFI_4', '6ZGD_4']:  # If pH 4...
                            for idx in range(0, 998):
                                Lx = superLambda[sim][rep][chain][idx]
                                Dx   = superData[sim][rep][chain][name][idx]

                                if Lx > 1 - protoC:  # If deprotonated...
                                    frameCount += 1
                                    if Dx < cutoff:  # If within cutoff...
                                        eventCount += 1

                    elif bin == 'proto_ph4':
                        for sim in ['4HFI_4', '6ZGD_4']:  # If pH 4...
                            for idx in range(0, 998):
                                Lx = superLambda[sim][rep][chain][idx]
                                Dx   = superData[sim][rep][chain][name][idx]

                                if Lx < protoC:  # If protonated...
                                    frameCount += 1
                                    if Dx < cutoff:  # If within cutoff...
                                        eventCount += 1

                    if frameCount == 0:
                        frac = 0
                    else:
                        frac = eventCount / float(frameCount)
                    superResult[bin][name].append(frac)

    #? MAKE THE BARPLOT

    # Initialize required data structures
    meanList = {}
    serrList = {}
    for key in bins:
        meanList[key] = []
        serrList[key] = []

    # Fill required data structures
    for bin in bins:
        for name in topList:
            meanList[bin].append(np.mean(superResult[bin][name]))
            serrList[bin].append(stderr(superResult[bin][name]))

    # print(meanList)  # debug

    # Some plotting stuff
    width = 0.2
    x     = np.arange(len(topList))
    fig   = plt.figure(figsize=(10, 6))
    ax    = fig.add_subplot()

    # deproto_ph7
    mean4 = meanList['deproto_ph7']
    serr4 = serrList['deproto_ph7']
    ax.bar(     x - width * 1.5, mean4, width, color='#8856a7', label='pH 7, deprotonated')
    ax.errorbar(x - width * 1.5, mean4, serr4, color='#8856a7', fmt='none', capsize=6, linewidth=2)

    # proto_ph7
    mean3 = meanList['proto_ph7']
    serr3 = serrList['proto_ph7']
    ax.bar(     x - width / 2.0, mean3, width, color='#9ebcda', label='pH 7, protonated')
    ax.errorbar(x - width / 2.0, mean3, serr3, color='#9ebcda', fmt='none', capsize=6, linewidth=2)

    # deproto_ph4
    mean2 = meanList['deproto_ph4']
    serr2 = serrList['deproto_ph4']
    ax.bar(     x + width / 2.0, mean2, width, color='#8856a7', label='pH 4, deprotonated', edgecolor='w', lw=1, hatch='//')
    ax.errorbar(x + width / 2.0, mean2, serr2, color='#8856a7', fmt='none', capsize=6, linewidth=2)

    # proto_ph4
    mean1 = meanList['proto_ph4']
    serr1 = serrList['proto_ph4']
    ax.bar(     x + width * 1.5, mean1, width, color='#9ebcda', label='pH 4, protonated', edgecolor='w', lw=1, hatch='\\\\')
    ax.errorbar(x + width * 1.5, mean1, serr1, color='#9ebcda', fmt='none', capsize=6, linewidth=2)

    ax.set_xticks(x, topList)
    ax.legend(loc=1, prop={'size': 12})

    # Do the title
    # u = MDAnalysis.Universe('../sims/4HFI_4/01/CA.pdb')
    # triplet = u.select_atoms('chainID A and resid {}'.format(target)).residues.resnames[0]
    # letter = triplet2letter(triplet)
    plt.title('Residue {}{} (structure mixing)'.format(letter, target))

    plt.ylim(0, 1.1)
    plt.ylabel('Contact occupancy')
    plt.tight_layout()
    plt.savefig('contacts/{}_strucmix.png'.format(target))
    plt.clf()
    plt.close()
