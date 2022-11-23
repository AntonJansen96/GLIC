#!/bin/python3

from science.utility import makeSuperDict

sims    = ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']
reps    = [1, 2, 3, 4]
chains  = ['A', 'B', 'C', 'D', 'E']
metrics = ['ecd_twist', 'ecd_spread', 'ecd_upper_spread', 'beta_expansion', 'm2_m1_dist', 'm2_radius', 'm1_kink', 'm1_kink_alt', 'nine_prime_dist', 'nine_prime_pore', 'minus_two_prime_dist', 'minus_two_prime_pore', 'c_loop']

superData = makeSuperDict([sims, reps, chains, metrics, []])

print(superData)
