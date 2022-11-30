#!/bin/python3

import json
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# Set global font size for figures.
matplotlib.rcParams.update({'font.size': 16})


def doPlot(basePath, name, label, sdev=True):
    # Load protein.json data file.
    data = json.load(open(f'{basePath}/{name}.json'))

    # Plot time-averaged pore radius.
    r = [i * 10 for i in data['pathwayProfile']['radiusMean']]
    z = [i * 10 for i in data['pathwayProfile']['s']]
    plt.plot(r, z, lw=2, label=label)

    # Plot associated standard deviation as a shaded region.
    if sdev:
        s = [i * 10 for i in data['pathwayProfile']['radiusSd']]
        A = np.array(z) + np.array(s)
        B = np.array(z) - np.array(s)
        plt.fill_between(r, A, B, facecolor='#000000', alpha=0.2)


for name in ['protein', 'backbone', 'c-alpha']:

    plt.figure(figsize=(5, 7))

    doPlot('6ZGD_static', name, '6ZGD', sdev=False)
    doPlot('4HFI_static', name, '4HFI', sdev=False)

    # for sim in ['4HFI_7', '4HFI_4']:
    for sim in ['6ZGD_7', '4HFI_7', '6ZGD_4', '4HFI_4']:
        doPlot(f'{sim}/comb', name, sim)

    plt.grid()
    plt.xlabel('Radius (Å)')
    plt.ylabel('Z-coordinate (Å)')
    plt.axis([0, 6, -20, 20])
    plt.legend(fontsize='small')
    plt.title(name)
    plt.tight_layout()
    plt.savefig(f'{name}.png')
    plt.clf()
