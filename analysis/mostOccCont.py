#!/bin/python3

import MDAnalysis
from science.utility import createIndexFile
from science.utility import gromacs
from science.utility import resname2triplet

# Input parameters
cutoff    = 0.40  # nm
pdb       = '../sims/4HFI_4/01/CA.pdb'
xtc       = '../sims/4HFI_4/01/MD_conv.xtc'

target    = 'resid 35 and chainid A and name OE1 OE2'
reference = 'protein'

# Create index file
createIndexFile(pdb, 'test.ndx', [reference, target])

# Run gmx select to obtain an index file containing the ATOM NUMBERS within
# the specified cutoff of the selection.
select = 'group 0 and within {} of group 1'.format(cutoff)
gromacs('select -f {} -n test.ndx -on out.ndx -select \'{}\''.format(xtc, select))

# Parse the out.ndx file into a list of frames where each frame is a list of atom numbers.
# We first parse the file into blocks, which are separated by an empty line.
# Then, within one block, we ignore the empty lines, and lines starting with '['.
# We add together any remaining lines (containing the atom numbers.)
frameList = []
for block in open('out.ndx').read().split('\n\n'):

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
with open('{}.txt'.format(target).replace(' ', '_'), 'w') as file:
    file.write('# Contact occupancy for \'{}\' (reference = \'{}\'), cutoff = {} nm\n'.format(target, reference, cutoff))

    for comb in bins:
        resLetter = resname2triplet(comb[2])
        file.write('{}{:<4d} {:2s} {:.4f}\n'.format(resLetter, comb[0], comb[1], bins[comb] / float(len(frameList))))

# add loop for different chains, sims, reps
# unroll this loop with multithreading
