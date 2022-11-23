#!/bin/python3

import json
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# Set global font size for figures.
matplotlib.rcParams.update({'font.size': 14})

sims = ['6ZGD_7', '6ZGD_4', '4HFI_7', '4HFI_4']
reps = [1]

plt.figure(figsize=(5, 7))

for sim in sims:
    for rep in reps:
        # Load protein.json data file.
        data = json.load(open('{}/{:02d}/protein.json'.format(sim, rep)))

        # Plot time-averaged pore radius.
        r = [i * 10 for i in data['pathwayProfile']['radiusMean']]
        z = [i * 10 for i in data['pathwayProfile']['s']]
        plt.plot(r, z, lw=2, label='{}_{}'.format(sim, rep))

        # Plot the associated standard deviation as a shaded region.
        s = [i * 10 for i in data['pathwayProfile']['radiusSd']]
        A = np.array(z) + np.array(s)
        B = np.array(z) - np.array(s)
        plt.fill_between(r, A, B, facecolor='#000000', alpha=0.2)

plt.xlabel('Radius (Å)')
plt.ylabel('Z-coordinate (Å)')
plt.axis([0, 6, -20, 20])
plt.legend(fontsize='small')
plt.title('built-in Sd (\'radiusSd\')')
plt.tight_layout()
plt.savefig('profile1.png')

#? DO THE STANDARD DEVATION/ERROR INSTEAD BY AVERAGING OVER FOUR REPLICAS:
