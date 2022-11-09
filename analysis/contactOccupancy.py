#!/bin/python3

import os
import MDAnalysis
import multiprocessing as mp

from science.utility import createIndexFile
from science.utility import gromacs
from science.utility import triplet2letter


def doContactOccupancy(pdb, xtc, target, ref='protein', outfile='', cutoff=0.40, idxA='select.ndx', idxB='output.ndx', outputDir='contacts'):
    """Performs contact occupancy analysis between and MDAnalysis-style selection
    called 'target' and a reference structure called 'ref'. Wraps gmx select.
    Final output is a list of most occupied contacts in descending order.

    Args:
        pdb (string): path to .pdb or .gro file.
        xtc (string): path to .trr or .xtc file.
        target (string): MDAnalysis style selection string.
        ref (str, optional): MDAnalysis style selection string. Defaults to 'protein'.
        outfile (str, optional): output file name. Defaults to 'target' string.
        cutoff (float, optional): contact cutoff distance in nm. Defaults to 0.40.
        idxA (str, optional): name of index file for 'ref' and 'target'. Defaults to 'select.ndx'.
        idxB (str, optional): name of index file for output of atoms within cutoff. Defaults to 'output.ndx'.
        outputDir (str, optional): directory for output. Defaults to 'contacts'.
    """

    # Prepare directory stuff
    if not os.path.exists(outputDir):
        os.mkdir(outputDir)

    idxA = outputDir + '/' + idxA
    idxB = outputDir + '/' + idxB

    # Create index file
    createIndexFile(pdb, idxA, [ref, target])

    # Run gmx select to obtain an index file containing the ATOM NUMBERS within
    # the specified cutoff of the selection.
    select = 'group 0 and within {} of group 1'.format(cutoff)
    gromacs('select -f {} -n {} -on {} -select \'{}\''.format(xtc, idxA, idxB, select))

    # Parse the idxB file into a list of frames where each frame is a list of atom numbers.
    # We first parse the file into blocks, which are separated by an empty line.
    # Then, within one block, we ignore the empty lines, and lines starting with '['.
    # We add together any remaining lines (containing the atom numbers.)
    frameList = []
    for block in open(idxB).read().split('\n\n'):

        totalString = ""
        for line in block.split('\n'):
            # Skip empty or header lines within the block
            if (not line) or (line[0] == '['):
                continue

            totalString += line

        # Add lines containing atom numbers to the frameList.
        frame = [int(val) for val in totalString.split()]
        frameList.append(frame)

    # Start binning per chain_residue.
    # Basically, we make a dictionary 'bins' with as keys tuples like (35, 'A', 'GLU')
    # and the value is the number of occurances of a contact in the frames.
    u = MDAnalysis.Universe(pdb)

    bins = {}
    for frame in frameList:

        contacts = []
        for atomNum in frame:

            resid   = u.atoms[atomNum - 1].resid
            chain   = u.atoms[atomNum - 1].chainID
            resname = u.atoms[atomNum - 1].resname
            # make the tuple and append to temporary list
            contacts.append((resid, chain, resname))
        # we use a set here as we may have contacts with multiple atoms in the same
        # residue, which DO NOT count double.
        for comb in set(contacts):
            if comb in bins:
                bins[comb] += 1     # if entry exists, simply increment.
            else:
                bins[comb] = 1      # if entry does not exist, initialize to 1.

    # Sort bins on descending number of contacts
    bins = dict(sorted(bins.items(), key=lambda item: item[1], reverse=True))

    # Write final output file.
    if outfile == '':
        file = open('{}/{}.txt'.format(outputDir, target).replace(' ', '_'), 'w')
    else:
        file = open('{}/{}'.format(outputDir, outfile), 'w')

    file.write('# Contact occupancy for \'{}\' (reference \'{}\')\n'.format(target, ref))
    file.write('# pdb = {}, xtc = {}, cutoff = {} nm\n'.format(pdb, xtc, cutoff))

    for comb in bins:
        resLetter = triplet2letter(comb[2])
        file.write('{}{:<4d} {:2s} {:.4f}\n'.format(resLetter, comb[0], comb[1], bins[comb] / float(len(frameList))))

    file.close()

    # Prevent too many GROMACS backups (error) from being made.
    os.system('rm {} {}'.format(idxA, idxB))


if __name__ == "__main__":

    # WRAPPER FOR MULTITHREADING
    def task(resid, sim, rep, chain):
        pdb     = '../sims/{}/{:02d}/CA.pdb'.format(sim, rep)
        xtc     = '../sims/{}/{:02d}/MD_conv.xtc'.format(sim, rep)
        target  = 'resid {} and chainid {} and name OE1 OE2 OD1 OD2 NE2 ND1'.format(resid, chain)
        outfile = '{}_{}_{}_{}.txt'.format(sim, rep, resid, chain)

        thread = mp.current_process().pid
        idxA   = 'select' + str(thread) + '.ndx'
        idxB   = 'output' + str(thread) + '.ndx'

        doContactOccupancy(pdb, xtc, target, outfile=outfile, idxA=idxA, idxB=idxB)

    # PREPARE ITERABLES
    items = []
    for resid in [13, 14, 26, 31, 32, 35, 49, 55, 67, 69, 75, 82, 86, 88, 91, 97, 104, 115, 122, 127, 136, 145, 147, 153, 154, 161, 163, 177, 178, 181, 185, 222, 235, 243, 272, 277, 282]:
        for sim in ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']:
            for rep in [1, 2, 3, 4]:
                for chain in ['A', 'B', 'C', 'D', 'E']:
                    items.append((resid, sim, rep, chain))

    # RUN MULTITHREADED
    pool = mp.Pool(processes=mp.cpu_count())
    pool.starmap(task, items, chunksize=1)
