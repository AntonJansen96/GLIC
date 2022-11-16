#!/bin/python3

import os
import multiprocessing as mp

from science.parsing import loadCol
from science.utility import gromacs
from science.utility import createIndexFile

outputDir = 'time'

# PREPARE DIRECTORY STUFF

if not os.path.exists(outputDir):
    os.mkdir(outputDir)

# WRAPPER(S) FOR MULTITHREADING


def task1(sim, rep, target, contact, chain1, chain2):
    """For the analysis of regular contacts.
    """
    # Unique threads need a unique index file name
    thread = mp.current_process().pid
    idxA   = outputDir + '/' + 'select' + str(thread) + '.ndx'

    # Create index file
    sel0 = 'resid {} and chainID {} and name OE1 OE2 OD1 OD2 NE2 ND1'.format(target, chain1)
    sel1 = 'resid {} and chainID {} and name NZ NE2 OD2 OE1 OD1 NE OG1 OH ND2 OE2 OT1 NH2 O NH1 NE1 N ND1 OT2 OG'.format(contact, chain2)
    createIndexFile('../sims/{}/{:02d}/CA.pdb'.format(sim, rep), idxA, [sel0, sel1])

    # Run gxm mindist to get the .xvg file for these selections
    p1 = "../sims/{}/{:02d}/CA.pdb".format(sim, rep)
    p2 = "../sims/{}/{:02d}/MD_conv.xtc".format(sim, rep)
    gromacs('mindist -s {} -f {} -n {} -od {}/{}_{}_{}{}_{}{}.xvg'.format(p1, p2, idxA, outputDir, sim, rep, target, chain1, contact, chain2), stdin=[0, 1])

    # Prevent too many GROMACS backups (error) from being made.
    os.system('rm {}'.format(idxA))


def task2(sim, rep, target, chain1):
    """For the analysis of Na+.
    """
    # Unique threads need a unique index file name
    thread = mp.current_process().pid
    idxA   = outputDir + '/' + str(thread) + '.ndx'

    # Create index file
    sel0 = 'resid {} and chainID {} and name OE1 OE2 OD1 OD2 NE2 ND1'.format(target, chain1)
    sel1 = 'resname NA'
    createIndexFile('../sims/{}/{:02d}/CA.pdb'.format(sim, rep), idxA, [sel0, sel1])

    # Run gxm mindist to get the .xvg file for these selections
    p1 = "../sims/{}/{:02d}/CA.pdb".format(sim, rep)
    p2 = "../sims/{}/{:02d}/MD_conv.xtc".format(sim, rep)

    gromacs('mindist -s {} -f {} -n {} -od {}/{}_{}_{}{}_NA.xvg'.format(p1, p2, idxA, outputDir, sim, rep, target, chain1), stdin=[0, 1])

    # Prevent too many GROMACS backups (error) from being made.
    os.system('rm {}'.format(idxA))


# PREPARE ITERABLES

items1 = []  # Holds tuples for regular contacts.
items2 = []  # Holds tuples specifically for Na+.

for sim in ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']:
    for rep in [1]:
        for target in [13, 14, 26, 31, 32, 35, 49, 55, 67, 69, 75, 82, 86, 88, 91, 97, 104, 115, 122, 127, 136, 145, 147, 153, 154, 161, 163, 177, 178, 181, 185, 222, 235, 243, 272, 277, 282]:

            # Load the top 6 contacts for said residue from the contacts/target.dat file
            topList = loadCol('contacts/{}.dat'.format(target), col=1, header=1)[:6]

            # If there is a Na+, remove it. If not, remove the last element.
            if 'Na+' in topList:
                topList.remove('Na+')
            else:
                topList = topList[:5]

            # Remove first letter from strings in topList
            topList = [val[1:] for val in topList]

            chain1 = ['A', 'B', 'C', 'D', 'E']

            for contact in topList:
                if contact[-1] == 'c':
                    chain2  = ['E', 'A', 'B', 'C', 'D']
                    contact = int(contact[:-1])

                elif contact[-1] == 'p':
                    chain2  = ['B', 'C', 'D', 'E', 'A']
                    contact = int(contact[:-1])

                else:
                    chain2  = ['A', 'B', 'C', 'D', 'E']
                    contact = int(contact)

                for idx in range(0, len(chain1)):
                    items1.append((sim, rep, target, contact, chain1[idx], chain2[idx]))

            for letter in chain1:
                items2.append((sim, rep, target, letter))

# RUN MULTITHREADED

pool = mp.Pool(processes=mp.cpu_count())
pool.starmap(task1, items1, chunksize=1)
pool.starmap(task2, items2, chunksize=1)
