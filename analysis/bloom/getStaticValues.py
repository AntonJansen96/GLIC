#!/bin/python3

import numpy as np

from science.parsing import loadCol
from science.utility import makeSuperDict

sims    = ['4HFI', '6ZGD']
chains  = ['A', 'B', 'C', 'D', 'E']
metrics = ['ecd_twist', 'ecd_spread', 'ecd_upper_spread', 'beta_expansion', 'm2_m1_dist', 'm2_radius', 'm1_kink', 'm1_kink_alt', 'nine_prime_dist', 'nine_prime_pore', 'minus_two_prime_dist', 'minus_two_prime_pore', 'c_loop']

superData = makeSuperDict([sims, metrics, []])

# Load the data.
for sim in sims:
    for chain in chains:
        fname = '{}_0_{}.txt'.format(sim, chain)
        for idx in range(0, len(metrics)):
            superData[sim][metrics[idx]] += loadCol(fname, idx + 1, header=0)

# Print result.
for sim in sims:
    print('\n{}'.format(sim))
    for metric in metrics:
        print('{}: {:.3f}'.format(metric, np.mean(superData[sim][metric])))
