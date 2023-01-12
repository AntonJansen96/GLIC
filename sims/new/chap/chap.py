#!/bin/python3

import json
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# Set global font size for figures.
matplotlib.rcParams.update({'font.size': 16})


# To fix a plotting bug.
def swap(*line_list):
    for lines in line_list:
        try:
            iter(lines)
        except:  # noqa
            lines = [lines]
        for line in lines:
            xdata, ydata = line.get_xdata(), line.get_ydata()
            line.set_xdata(ydata)
            line.set_ydata(xdata)
            line.axes.autoscale_view()


def doPlot(basePath, name, label, color, style, sdev):
    # Load protein.json data file.
    data = json.load(open(f'{basePath}/{name}.json'))

    # Plot time-averaged pore radius.
    x = [i * 10 for i in data['pathwayProfile']['s']]
    y = [i * 10 for i in data['pathwayProfile']['radiusMean']]

    line1 = plt.plot(x, y, lw=1.5, label=label, color=color, linestyle=style)
    swap(line1)

    # Plot associated standard deviation as a shaded region.
    if sdev:
        errors = [i * 10 for i in data['pathwayProfile']['radiusSd']]

        A = np.array(y) + np.array(errors)
        B = np.array(y) - np.array(errors)

        plt.fill_betweenx(x, A, B, facecolor=color, alpha=0.1)


################################################################################

# PROTEIN, FIGURE 1

plt.figure(figsize=(5, 7))

for type in ['protein', 'backbone']:

    doPlot('4HFI_static', type, '4HFI crystal', color='black',   style='-',  sdev=True)
    doPlot('4HFI_4',      type, '4HFI_4',       color='#1f77b4', style='--', sdev=True)
    doPlot('4HFI_4_new',  type, '4HFI_4 New',   color='#1f77b4', style='-',  sdev=True)
    doPlot('4HFI_7',      type, '4HFI_7',       color='#ff7f0e', style='--', sdev=True)
    doPlot('4HFI_7_new',  type, '4HFI_7 New',   color='#ff7f0e', style='-',  sdev=True)

    if type == 'protein':
        plt.axis([0, 9, -20, 50])
    elif type == 'backbone':
        plt.axis([0, 12, -20, 50])

    plt.grid()
    plt.xlabel('Radius (Å)')
    plt.ylabel('Z-coordinate (Å)')
    plt.legend(fontsize='small')
    plt.title(type)
    plt.tight_layout()
    plt.savefig(f'{type}1.png')
    plt.clf()

for type in ['protein', 'backbone']:

    doPlot('6ZGD_static', type, '6ZGD cryoEM',  color='black',   style='-',  sdev=False)
    doPlot('6ZGD_7',      type, '6ZGD_7',       color='#8856a7', style='-',  sdev=True)
    doPlot('6ZGD_4',      type, '6ZGD_4',       color='#9ebcda', style='-',  sdev=True)
    doPlot('4HFI_static', type, '4HFI crystal', color='black',   style='--', sdev=False)
    doPlot('4HFI_7',      type, '4HFI_7',       color='#8856a7', style='--', sdev=True)
    doPlot('4HFI_4',      type, '4HFI_4',       color='#9ebcda', style='--', sdev=True)

    if type == 'protein':
        plt.axis([0, 9, -20, 50])
    elif type == 'backbone':
        plt.axis([0, 12, -20, 50])

    plt.grid()
    plt.xlabel('Radius (Å)')
    plt.ylabel('Z-coordinate (Å)')
    plt.legend(fontsize='small')
    plt.title(type)
    plt.tight_layout()
    plt.savefig(f'{type}2.png')
    plt.clf()
