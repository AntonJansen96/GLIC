#!/usr/bin/env python3

import pickle
import numpy as np

sims = ['6ZGD_7', '6ZGD_4', '4HFI_7', '4HFI_4']
letters = ['D', 'E', 'E', 'D', 'D', 'E', 'D', 'D', 'E', 'E', 'E', 'E', 'D', 'D', 'D', 'D', 'E', 'D', 'D', 'H', 'D', 'D', 'E', 'D', 'D', 'D', 'E', 'E', 'D', 'E', 'D', 'E', 'H', 'E', 'E', 'H', 'E']
residues = [13, 14, 26, 31, 32, 35, 49, 55, 67, 69, 75, 82, 86, 88, 91, 97, 104, 115, 122, 127, 136, 145, 147, 153, 154, 161, 163, 177, 178, 181, 185, 222, 235, 243, 272, 277, 282]

protoMean = pickle.load(open('paperBarPlots.obj', 'rb'))

with open('protonation.txt', 'w+') as file:
    for sim in sims:
        for idx in range(0, len(residues)):
            protoMean[sim][residues[idx]] = np.mean(protoMean[sim][residues[idx]])
            file.write(f"{sim} {letters[idx]}{residues[idx]} {protoMean[sim][residues[idx]]:.3f}\n")

################################################################################

cat1 = []
cat2 = []
cat3 = []

cutoff = 0.20  # 20%.

low  = {'D': 0.41, 'E': 0.56, 'H': 0.93}
high = {'D': 0.03, 'E': 0.06, 'H': 0.38}

for idx in range(0, len(residues)):
    for sim in sims:

        if sim in ['4HFI_4', '6ZGD_4']:
            theo = low
        else:
            theo = high

        if abs(protoMean[sim][residues[idx]] - theo[letters[idx]]) > cutoff:
            cat2.append(f"{letters[idx]}{residues[idx]}")
            break

    if f"{letters[idx]}{residues[idx]}" not in cat2:
        cat1.append(f"{letters[idx]}{residues[idx]}")

print('category 1', cat1)
print('category 2', cat2)

# See https://docs.google.com/spreadsheets/d/1sS74Zz6NqjqPqtH6sFTLH-WD5Ek0QppPBB0BPUbGjSg/edit#gid=41793343
# Cutoff 25%

# CAT 1 (11): no pH or state dependence
# D13, E14, E67, E69, E104, H127, E147, D154, D161, D178, E222

# CAT 2 (20): pH-dependence but not state dependence
# D31, D32, E35, D49, D55, E75, D86, D88, D91, D115, D136, D145, D153, E163, D185, H235, E243, E272, H277, E282

# CAT 3 (6): state dependence (and additionally pH-dependence)
# E26, E82, D97, D122, E177, E181
