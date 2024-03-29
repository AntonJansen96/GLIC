#!/bin/python3

import matplotlib
import matplotlib.pyplot as plt
import os
import numpy as np
import ternary
import pickle
import sys

from science.parsing import loadxvg
from science.parsing import Structure
from science.utility import Stopwatch
from science.cphmd import protonation

from data import biophys, nury2010, fritsch2011, lev2017, nemecz2017

# Set global font size for figures.
matplotlib.rcParams.update({'font.size': 14})


class TwoState:
    """Holds data for a two-state (ASP or GLU) titratable group."""
    def __init__(self, idx, resname, resid, chain, t, x):
        self.d_idx      = idx
        self.d_fname    = 'cphmd-coord-{}.xvg'.format(idx)
        self.d_resname  = resname
        self.d_resid    = resid
        self.d_chain    = chain
        self.d_t        = t
        self.d_x        = x


class MultiState:
    """Holds data for a multistate-state (HIS) titratable group."""
    def __init__(self, idx, resname, resid, chain, t1, t2, t3, x1, x2, x3):
        self.d_idx      = [idx, idx + 1, idx + 2]
        self.d_fname    = ['cphmd-coord-{}.xvg'.format(idx), 'cphmd-coord-{}.xvg'.format(idx + 1), 'cphmd-coord-{}.xvg'.format(idx + 2)]
        self.d_resname  = resname
        self.d_resid    = resid
        self.d_chain    = chain
        self.d_t        = [t1, t2, t3]
        self.d_x        = [x1, x2, x3]


class Buffer:
    """Holds data for the buffer group."""
    def __init__(self, idx, t, x, count):
        self.d_idx     = idx
        self.d_fname   = 'cphmd-coord-{}.xvg'.format(idx)
        self.d_resname = 'BUF'
        self.d_t       = t
        self.d_x       = x
        self.d_count   = count


class Replica:
    """Holds data for one replica."""
    def __init__(self, name, replicaID, twoStateList, multiStateList, buffer):
        self.d_name           = name
        self.d_replicaID      = replicaID
        self.d_twoStateList   = twoStateList
        self.d_multiStateList = multiStateList
        self.d_buffer         = buffer


class ReplicaSet:
    """Holds data for one set of GLIC replicas (e.g. 4HFI_4)."""
    def __init__(self, name, replicaSet):
        self.d_name        = name
        self.d_replica     = replicaSet


class GLICSims:
    """Holds the data for all four GLIC replica sets"""
    def __init__(self, folders, replicas, dt, b):
        os.chdir('../sims')

        self.d_replicaSet = []
        for ii in folders:
            os.chdir(ii)

            replicaList = []
            for jj in replicas:
                os.chdir(jj)

                idx = 1
                foundBUF = False
                twoStateList   = []
                multiStateList = []
                for residue in Structure('CA.pdb', 2).d_residues:

                    if residue.d_resname in ['ASPT', 'GLUT']:
                        print('{} : {} Loading {}-{} in chain {}...'.format(ii, jj, residue.d_resname, residue.d_resid, residue.d_chain), end='\r')

                        xvgdata = loadxvg('cphmd-coord-{}.xvg'.format(idx), dt=dt, b=b)

                        twoStateList.append(TwoState(
                            idx,
                            residue.d_resname,
                            residue.d_resid,
                            residue.d_chain,
                            xvgdata[0],  # time
                            xvgdata[1]   # coordinates
                        ))

                        idx += 1

                    elif residue.d_resname == 'HSPT':
                        print('{} : {} Loading {}-{} in chain {}...'.format(ii, jj, residue.d_resname, residue.d_resid, residue.d_chain), end='\r')

                        xvgdata1 = loadxvg('cphmd-coord-{}.xvg'.format(idx),     dt=dt, b=b)
                        xvgdata2 = loadxvg('cphmd-coord-{}.xvg'.format(idx + 1), dt=dt, b=b)
                        xvgdata3 = loadxvg('cphmd-coord-{}.xvg'.format(idx + 2), dt=dt, b=b)

                        multiStateList.append(MultiState(
                            idx,
                            residue.d_resname,
                            residue.d_resid,
                            residue.d_chain,
                            xvgdata1[0],  # file 1 time
                            xvgdata2[0],  # file 2 time
                            xvgdata3[0],  # file 3 time
                            xvgdata1[1],  # file 1 coordinates
                            xvgdata2[1],  # file 2 coordinates
                            xvgdata3[1])  # file 3 coordinates
                        )

                        idx += 3

                    elif residue.d_resname == 'BUF' and not foundBUF:
                        xvgdata  = loadxvg('cphmd-coord-{}.xvg'.format(idx), dt=dt)
                        buffer   = Buffer(idx, xvgdata[0], xvgdata[1], count=185)
                        foundBUF = True

                replicaList.append(Replica(ii, int(jj), twoStateList, multiStateList, buffer))

                os.chdir('..')

            self.d_replicaSet.append(ReplicaSet(ii, replicaList))

            os.chdir('..')

        os.chdir('../analysis')

        # DIRECTORY STRUCTURE
        if not os.path.isdir('lambdaplots'):
            os.mkdir('lambdaplots')

    # FUNCTION FOR MOVING DEPROTONATION
    def movingDeprotonation(self, tList, xList, window):
        Sum = sum(xList[0:window])  # Initialize

        t = tList[window:]
        x = len(range(window, len(xList))) * [0]

        for ii in range(window, len(xList)):
            x[ii - window] = Sum / float(window)
            Sum -= xList[ii - window]
            Sum += xList[ii]

        return t, x

    def histograms(self, b):
        totalResidues = len(self.d_replicaSet[0].d_replica[0].d_twoStateList)
        chains = 5
        residuesPerChain = int(totalResidues / chains)

        # Outer most loop is over the ReplicaSets (4HFI_4, 4HFI_7, etc.):
        for ii in range(0, len(self.d_replicaSet)):
            # Second loop is over the titratable residues:
            for jj in range(0, residuesPerChain):
                # Set valuesList and protoFracList to zero:
                valuesList = []
                protoFracList = []
                # Third loop is over the four replicas:
                for kk in range(0, len(self.d_replicaSet[ii].d_replica)):  # 4 replicas...
                    # And fourth loop is over the five chains:
                    for ll in range(0, chains):  # ...x5 chains = 20 samples

                        # GET THE DATA
                        x = self.d_replicaSet[ii].d_replica[kk].d_twoStateList[jj + residuesPerChain * ll].d_x
                        x = x[b:]  # Only analyse starting from the b=frame.
                        x = [1.0 - val for val in x]  # Mirror in vertical x=0.5 axis

                        # GET PROTONATION FRACTION
                        protoFracList.append(protonation(x))

                        # GET HISTOGRAM VALUES
                        values, bins = np.histogram(x, density=True, bins=200, range=(-0.1, 1.1))
                        valuesList.append(values)

                # COMPUTE MEAN AND STANDARD ERROR
                meanList  = len(values) * [0]  # 200, to hold mean for each bin
                errorList = len(values) * [0]  # 200, to hold erro for each bin

                for kk in range(0, len(values)):  # 200

                    # Create list of 20 values
                    temp = [0] * len(valuesList)  # 4*5=20
                    for ll in range(0, len(valuesList)):  # 4*5=20
                        temp[ll]  = valuesList[ll][kk]

                    meanList[kk]  = np.mean(temp)
                    errorList[kk] = np.std(temp)

                # PLOT MEAN AND SHADED REGION (ERROR)
                A = []
                B = []
                for kk in range(0, len(meanList)):
                    A.append(meanList[kk] + errorList[kk])
                    B.append(meanList[kk] - errorList[kk])

                plt.figure(figsize=(8, 6))
                plt.plot(bins[1:], meanList)
                plt.fill_between(bins[1:], A, B, alpha=0.4, color='#1f77b4')

                # MAKE PLOT MORE NICE
                plt.text(x=1, y=17.7, s='Proto\n$q = 0$',    ha='center', fontsize=12)
                plt.text(x=0, y=17.6, s='Deproto\n$q = -1$', ha='center', fontsize=12)
                plt.title(self.d_replicaSet[ii].d_name, fontsize=18)
                plt.axis([-0.1, 1.1, -0.1, 17.5])
                plt.xlabel(r"$\lambda$-coordinate")
                plt.xticks(ticks=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0], labels=[1.0, 0.8, 0.6, 0.4, 0.2, 0.0])  # because we mirror in vertical x=0.5 axis
                plt.grid()

                group = self.d_replicaSet[0].d_replica[0].d_twoStateList[jj]

                # PRINT LIST OF PROTONATIONS TO TEXT FILE
                with open('lambdaplots/{}_protonation.txt'.format(self.d_replicaSet[ii].d_name), 'a') as file:
                    file.write('{:03d}-{} {:.2f}\n'.format(group.d_resid, group.d_resname, np.mean(protoFracList)))

                # ADD EXPERIMENTAL DATA FOR PH=4 CASE
                if float(self.d_replicaSet[ii].d_name[5:]) == 4.0:
                    plt.vlines(x=biophys["{0}-{1}".format(group.d_resname, group.d_resid)], ymin=0, ymax=10, color='r', linewidth=4.0, label="biophysics.se/Prevost2012")
                    plt.vlines(x=nury2010["{0}-{1}".format(group.d_resname, group.d_resid)], ymin=0, ymax=8, color='g', linewidth=4.0, label="Nury2010/Cheng2010/Calimet2013")
                    plt.vlines(x=fritsch2011["{0}-{1}".format(group.d_resname, group.d_resid)], ymin=0, ymax=6, color='b', linewidth=4.0, label="Fritsch2011")
                    plt.vlines(x=lev2017["{0}-{1}".format(group.d_resname, group.d_resid)], ymin=0, ymax=4, color='c', linewidth=4.0, label="Lev2017")
                    plt.vlines(x=nemecz2017["{0}-{1}".format(group.d_resname, group.d_resid)], ymin=0, ymax=2, color='m', linewidth=4.0, label="Nemecz2017/Hu2018")
                    plt.legend(loc='upper center')

                # SAVE AND CLEAR
                plt.tight_layout()
                plt.savefig('lambdaplots/{}_{:03d}-{}.png'.format(self.d_replicaSet[ii].d_name, group.d_resid, group.d_resname))
                plt.clf()

    def convergence(self):
        totalResidues = len(self.d_replicaSet[0].d_replica[0].d_twoStateList)
        chains = 5
        residuesPerChain = int(totalResidues / chains)

        # Outer most loop is over the ReplicaSets (4HFI_4, 4HFI_7, etc.):
        for ii in range(0, len(self.d_replicaSet)):
            # Second loop is over the titratable residues:
            for jj in range(0, residuesPerChain):
                # Set valuelist to zero.
                valuesList = []
                # Third loop is over the four replicas:
                for kk in range(0, len(self.d_replicaSet)):  # 4 replicas...
                    # And fourth loop is over the five chains:
                    for ll in range(0, chains):  # ...x5 chains = 20 samples
                        # GET THE DATA
                        t = self.d_replicaSet[ii].d_replica[kk].d_twoStateList[jj + residuesPerChain * ll].d_t
                        x = self.d_replicaSet[ii].d_replica[kk].d_twoStateList[jj + residuesPerChain * ll].d_x
                        x = [1.0 - val for val in x]  # Mirror in vertical x=0.5 axis
                        valuesList.append(x)

                # COMPUTE MEAN AND STANDARD ERROR
                meanList  = []
                errorList = []

                for kk in range(0, len(x)):  # 200

                    # Create list of 20 values
                    temp = [0] * len(valuesList)  # 4*5=20
                    for ll in range(0, len(valuesList)):  # 4*5=20
                        temp[ll]  = valuesList[ll][kk]

                    meanList.append(np.mean(temp))
                    errorList.append(np.std(temp))

                # PLOT MEAN AND SHADED REGION (ERROR)
                A = []
                B = []
                for kk in range(0, len(meanList)):
                    A.append(meanList[kk] + errorList[kk])
                    B.append(meanList[kk] - errorList[kk])

                t = [val * 0.001 for val in t]

                plt.figure(figsize=(8, 6))
                plt.plot(t, meanList)
                plt.fill_between(t, A, B, alpha=0.4, color='#1f77b4')

                # # PLOT
                # a, b = self.movingDeprotonation(t, x, window)
                # plt.plot(a, b)

                # MAKE PLOT MORE NICE
                plt.title(self.d_replicaSet[ii].d_name, fontsize=18)
                plt.ylim(-0.1, 1.1)
                plt.xlabel("Time (ns)")
                plt.ylabel("Protonation")
                plt.grid()

                group = self.d_replicaSet[0].d_replica[0].d_twoStateList[jj]

                # SAVE AND CLEAR
                plt.tight_layout()
                plt.savefig('lambdaplots/{}_{:03d}-{}_conv.png'.format(self.d_replicaSet[ii].d_name, group.d_resid, group.d_resname))
                plt.clf()

    def convergence_old(self, window):
        totalResidues = len(self.d_replicaSet[0].d_replica[0].d_twoStateList)
        chains = 5
        residuesPerChain = int(totalResidues / chains)

        # Outer most loop is over the ReplicaSets (4HFI_4, 4HFI_7, etc.):
        for ii in range(0, len(self.d_replicaSet)):
            # Second loop is over the titratable residues:
            for jj in range(0, residuesPerChain):
                # Third loop is over the four replicas:
                for kk in range(0, len(self.d_replicaSet[ii].d_replica)):  # 4 replicas...
                    # And fourth loop is over the five chains:
                    chainName = ['A', 'B', 'C', 'D', 'E']
                    for ll in range(0, chains):  # ...x5 chains = 20 samples

                        # GET THE DATA
                        t = self.d_replicaSet[ii].d_replica[kk].d_twoStateList[jj + residuesPerChain * ll].d_t
                        x = self.d_replicaSet[ii].d_replica[kk].d_twoStateList[jj + residuesPerChain * ll].d_x
                        x = [1.0 - val for val in x]  # Mirror in vertical x=0.5 axis

                        # PLOT
                        a, b = self.movingDeprotonation(t, x, window)
                        a = [val * 0.001 for val in a]
                        plt.plot(a, b, label=chainName[ll])

                # MAKE PLOT MORE NICE
                plt.title(self.d_replicaSet[ii].d_name, fontsize=18)
                plt.ylim(-0.1, 1.1)
                plt.xlabel("Time (ns)")
                plt.ylabel("Protonation running average")
                plt.grid()
                # plt.legend()

                group = self.d_replicaSet[0].d_replica[0].d_twoStateList[jj]

                # SAVE AND CLEAR
                plt.tight_layout()
                plt.savefig('lambdaplots/{}_{:03d}-{}_conv.png'.format(self.d_replicaSet[ii].d_name, group.d_resid, group.d_resname))
                plt.clf()

    def histidine(self, b):
        totalResidues = len(self.d_replicaSet[0].d_replica[0].d_multiStateList)
        chains = 5
        residuesPerChain = int(totalResidues / chains)

        # Outer most loop is over the ReplicaSets (4HFI_4, 4HFI_7, etc.):
        for ii in range(0, len(self.d_replicaSet)):
            # Second loop is over the titratable residues (01, 02, etc.):
            for jj in range(0, residuesPerChain):
                # Start figure here
                plt.figure(figsize=(8, 6))
                # Third loop is over the three lambda groups:
                for kk in [0, 1, 2]:
                    # Set valuesList to zero.
                    valuesList = []
                    # Fourth loop is over the four replicas:
                    for ll in range(0, len(self.d_replicaSet[ii].d_replica)):  # 4 replicas...
                        # And fifth loop is over the five chains:
                        for mm in range(0, chains):  # ...x5 chains = 20 samples

                            # GET THE DATA
                            x = self.d_replicaSet[ii].d_replica[kk].d_multiStateList[jj + residuesPerChain * ll].d_x[kk]
                            x = x[b:]  # Only analyse starting from the b=frame.
                            x = [1.0 - val for val in x]  # Mirror in vertical x=0.5 axis

                            # GET HISTOGRAM VALUES, BINS
                            values, bins = np.histogram(x, density=True, bins=200, range=(-0.1, 1.1))
                            valuesList.append(values)

                    # COMPUTE MEAN AND STANDARD ERROR
                    meanList  = len(values) * [0]  # 200, to hold mean for each bin
                    errorList = len(values) * [0]  # 200, to hold error for each bin

                    for ll in range(0, len(values)):  # 200

                        # Create list of 20 values
                        temp = [0] * len(valuesList)  # 4*5=20
                        for mm in range(0, len(valuesList)):  # 4*5=20
                            temp[mm] = valuesList[mm][ll]

                        meanList[ll] = np.mean(temp)
                        errorList[ll] = np.std(temp)

                    # PLOT MEAN AND SHADED REGION (ERROR)
                    A = []
                    B = []
                    for ll in range(0, len(meanList)):
                        A.append(meanList[ll] + errorList[ll])
                        B.append(meanList[ll] - errorList[ll])

                    description = ['state 1 (double proto)', 'state 2 (anti)', 'state 3 (syn)']
                    color = ['#1f77b4', '#ff7f0e', '#2ca02c']

                    plt.plot(bins[1:], meanList, color=color[kk], label=description[kk])
                    plt.fill_between(bins[1:], A, B, alpha=0.4, color=color[kk])

                # MAKE PLOT MORE NICE
                plt.title(self.d_replicaSet[ii].d_name, fontsize=18)
                plt.axis([-0.1, 1.1, -0.1, 17.5])
                plt.xlabel(r"$\lambda$-coordinate")
                plt.xticks(ticks=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0], labels=[1.0, 0.8, 0.6, 0.4, 0.2, 0.0])  # because we mirror in vertical x=0.5 axis
                plt.grid()
                plt.legend(loc='upper center')

                group = self.d_replicaSet[0].d_replica[0].d_multiStateList[jj]

                # SAVE AND CLEAR
                plt.tight_layout()
                plt.savefig('lambdaplots/{}_{:03d}-{}.png'.format(self.d_replicaSet[ii].d_name, group.d_resid, group.d_resname))
                plt.clf()

    def hisheatmap(self, b):
        totalResidues = len(self.d_replicaSet[0].d_replica[0].d_multiStateList)
        chains = 5
        residuesPerChain = int(totalResidues / chains)

        # Outer most loop is over the ReplicaSets (4HFI_4, 4HFI_7, etc.):
        for ii in range(0, len(self.d_replicaSet)):
            # Second loop is over the titratable residues (01, 02, etc.):
            for jj in range(0, residuesPerChain):
                # Holds lambda1 and lambda2
                lambdaList = [[], []]
                # Third loop is over the three lambda groups:
                for kk in [0, 1]:
                    # Fourth loop is over the four replicas:
                    for ll in range(0, len(self.d_replicaSet[ii].d_replica)):  # 4 replicas...
                        # And fifth loop is over the five chains:
                        for mm in range(0, chains):  # ...x5 chains = 20 samples
                            # GET THE DATA
                            x = self.d_replicaSet[ii].d_replica[kk].d_multiStateList[jj + residuesPerChain * ll].d_x[kk]
                            x = x[b:]  # Only analyse starting from the b=frame.
                            lambdaList[kk] += x

                # MAKE HISTOGRAM
                Nbins = 11
                lambda1 = lambdaList[0]
                lambda2 = lambdaList[1]

                # Fix slightly uneven length problem
                if len(lambda1) > len(lambda2):
                    lambda1 = lambda1[0:len(lambda2)]
                elif len(lambda1) < len(lambda2):
                    lambda2 = lambda2[0:len(lambda1)]

                H, _ = np.histogramdd((lambda1, lambda2), bins=(Nbins, Nbins), range=((-0.1, 1.1), (-0.1, 1.1)), density=True)

                # code from Pavel I don't understand but is necessary
                interp_dict = dict()
                binx = np.linspace(0, 1, Nbins)
                for i, _ in enumerate(binx):
                    for j, _ in enumerate(binx):
                        interp_dict[(i, j)] = H[i, j]

                # TRIANGLE PLOTTING USING TERNARY
                fig, tax = ternary.figure(scale=Nbins - 1)
                tax.heatmap(interp_dict, cmap="OrRd", colorbar=False, vmin=0, vmax=20)

                # MAKE PLOT MORE NICE
                tax.boundary(linewidth=1.0)
                tax.set_title('{}\n'.format(self.d_replicaSet[ii].d_name), fontsize=13)
                tax.get_axes().axis('off')
                tax.clear_matplotlib_ticks()
                tax.right_corner_label("double", fontsize=11, offset=0.15)
                tax.top_corner_label("anti", fontsize=11, offset=0.15)
                tax.left_corner_label("   syn", fontsize=11, offset=0.15)
                fig.gca().set_aspect('equal')
                group = self.d_replicaSet[0].d_replica[0].d_multiStateList[jj]

                # SAVE, TRIM, AND RESIZE
                fname = 'lambdaplots/heat_{}_{:03d}-{}.png'.format(self.d_replicaSet[ii].d_name, group.d_resid, group.d_resname)
                tax.savefig(fname)
                os.system('convert {} -trim {}'.format(fname, fname))  # trim all whitespace
                os.system('convert -bordercolor white -border 6 {} {}'.format(fname, fname))  # add 6 pixels back on all sides
                os.system('convert {} -resize 73% {}'.format(fname, fname))  # resize so it fits with other HSPT histograms

    def doFinalPlots(self):
        os.chdir('lambdaplots')
        for res in ['127-HSPT', '235-HSPT', '277-HSPT']:
            # CREATE *OLD* FINAL HISTOGRAMS FOR HISTIDINE
            os.system('convert 6ZGD_7_{}.png 4HFI_7_{}.png +append temp1.png'.format(res, res))
            os.system('convert 6ZGD_4_{}.png 4HFI_4_{}.png +append temp2.png'.format(res, res))
            os.system('convert temp1.png temp2.png -append hist_{}.png'.format(res))
            # CREATE *NEW* FINAL HISTOGRAMS FOR HISTIDINE
            os.system('convert heat_6ZGD_7_{}.png heat_4HFI_7_{}.png +append temp1.png'.format(res, res))
            os.system('convert heat_6ZGD_4_{}.png heat_4HFI_4_{}.png +append temp2.png'.format(res, res))
            os.system('convert temp1.png temp2.png -append heat_{}.png'.format(res))
            # MERGE BOTH TO CREATE FINAL PLOTS
            os.system('convert hist_{}.png heat_{}.png +append final_{}.png'.format(res, res, res))

        for res in ['013-ASPT', '014-GLUT', '026-GLUT', '031-ASPT', '032-ASPT', '035-GLUT', '049-ASPT', '055-ASPT', '067-GLUT', '069-GLUT', '075-GLUT', '082-GLUT', '086-ASPT', '088-ASPT', '091-ASPT', '097-ASPT', '104-GLUT', '115-ASPT', '122-ASPT', '136-ASPT', '145-ASPT', '147-GLUT', '153-ASPT', '154-ASPT', '161-ASPT', '163-GLUT', '177-GLUT', '178-ASPT', '181-GLUT', '185-ASPT', '222-GLUT', '243-GLUT', '272-GLUT', '282-GLUT']:
            # CREATE FINAL HISTOGRAMS FOR ASPARTIC AND GLUTAMIC ACID
            os.system('convert 6ZGD_7_{}.png 4HFI_7_{}.png +append temp1.png'.format(res, res))
            os.system('convert 6ZGD_4_{}.png 4HFI_4_{}.png +append temp2.png'.format(res, res))
            os.system('convert temp1.png temp2.png -append hist_{}.png'.format(res))
            # CREATE FINAL CONVERGENCE PLOTS FOR ASPARTIC AND GLUTAMIC ACID
            os.system('convert 6ZGD_7_{}_conv.png 4HFI_7_{}_conv.png +append temp1.png'.format(res, res))
            os.system('convert 6ZGD_4_{}_conv.png 4HFI_4_{}_conv.png +append temp2.png'.format(res, res))
            os.system('convert temp1.png temp2.png -append conv_{}.png'.format(res))
            # MERGE BOTH TO CREATE FINAL PLOTS
            os.system('convert hist_{}.png conv_{}.png +append final_{}.png'.format(res, res, res))


if __name__ == "__main__":

    if sys.argv[1] == 'real':

        dt = 20         # Only use 1/20 frames.
        b  = 0          # Do not ignore any frames, load entire trajectory.
        # b  = 300000   # Ignore first 300'000 steps = 300 ns.

    elif sys.argv[1] == 'test':

        os.chdir('test')
        dt = 1
        b  = 2000

    else:
        print("arguments are \'real\' or \'test\'")
        exit()

    # PICKLE PART

    if not os.path.isfile('data.pickle'):

        print('No data.pickle detected, will construct and dump GLIC opbject...')

        t = Stopwatch('Constructing GLIC object')
        GLIC = GLICSims(['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7'], ['01', '02', '03', '04'], dt=dt, b=b)
        t.time()

        t = Stopwatch('Dumping GLIC object')
        with open('data.pickle', 'wb') as handle:
            pickle.dump(GLIC, handle, protocol=pickle.HIGHEST_PROTOCOL)
        t.time()

    else:
        print('Loading GLIC object from pickle data to save time...')

        t = Stopwatch('Loading GLIC object from data')
        with open('data.pickle', 'rb') as handle:
            GLIC = pickle.load(handle)
        t.time()

    # ACTUAL ANALYSIS CODE

    t = Stopwatch('Making plots')

    GLIC.histograms(b=15000)
    GLIC.convergence()
    GLIC.histidine(b=15000)
    GLIC.hisheatmap(b=15000)
    GLIC.doFinalPlots()

    t.time()
