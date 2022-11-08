#!/bin/python3

import os
import MDAnalysis
from science.utility import createIndexFile
from science.utility import gromacs
from science.utility import resname2triplet


def doContactOccupancy(pdb, xtc, target, ref='protein', outfile='', cutoff=0.40, idxA='select.ndx', idxB='output.ndx', outputDir='contacts'):
    """Performs contact occupancy analysis between and MDAnalysis-style selection
    called 'target' and a reference structure called 'ref'. Wraps gmx select.
    Final output is a list of most occupied contacts in descending order.

    Args:
        pdb (string): path to .pdb or .gro file.
        xtc (string): path to .trr or .xtc file.
        target (string): MDAnalysis style selection string.
        ref (str, optional): MDAnalysis style selection string. Defaults to 'protein'.
        outfile (str, optional): output file name. Defaults to target string.
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
        resLetter = resname2triplet(comb[2])
        file.write('{}{:<4d} {:2s} {:.4f}\n'.format(resLetter, comb[0], comb[1], bins[comb] / float(len(frameList))))

    file.close()


if __name__ == "__main__":

    pdb       = '../sims/4HFI_4/01/CA.pdb'
    xtc       = '../sims/4HFI_4/01/MD_conv.xtc'
    target    = 'resid 35 and chainid A and name OE1 OE2'

    doContactOccupancy(pdb, xtc, target)
