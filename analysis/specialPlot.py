#!/bin/python3

import matplotlib
import copy

from science.parsing import loadCol
from science.parsing import loadxvg
from science.cphmd import getLambdaFileIndices


# Set global font size for figures.
matplotlib.rcParams.update({'font.size': 14})

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


sims    = ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']
reps    = [1]
chains  = ['A', 'B', 'C', 'D', 'E']
targets = [35]

for target in targets:

    lambdaIndices = getLambdaFileIndices('../sims/4HFI_4/01/CA.pdb', target)

    #? CONSTRUCT THE SUPERDATA
    # (dict of sims of reps of chains of topResidues of (mindist) lists)

    # Load the top 6 contacts for said residue from the contacts/target.dat file
    topList = loadCol('contacts/{}.dat'.format(target), col=1, header=1)[:6]
    # If there is a Na+, remove it. If not, remove the last element.
    # Then, include Na+ by default at the end.
    if 'Na+' in topList:
        topList.remove('Na+')
    else:
        topList = topList[:5]
    topList.append('Na+')

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

    print(superData)

    #? CONSTRUCT THE SUPERLAMBDA
    # (dict of sims of reps of chains of lists of coordinates)

    A = {}
    for key in chains:
        A[key] = []
    B = {}
    for key in reps:
        B[key] = copy.deepcopy(A)
    superLambda = {}
    for key in sims:
        superLambda[key] = copy.deepcopy(B)

    #? GATHER ALL THE DATA INTO SUPERDATA AND SUPERLAMBDA
    for sim in sims:
        for rep in reps:
            for chain in chains:

                #! Fill superLambda
                fname = '../sims/{}/{:02d}/cphmd-coord-{}.xvg'.format(sim, rep, lambdaIndices[chains.index(chain)])
                superLambda[sim][rep][chain] = loadxvg(fname, dt=1000)[1]
                print(fname, chain)  # debug

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

                    print(fname)  # debug

    #? DO THE BINNING ETC
    