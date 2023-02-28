#!/usr/bin/env python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from science.cphmd import theoreticalProtonation
from science.utility import makeSuperDict

matplotlib.rcParams.update({'font.size': 20})

sims     = ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']
resnames = ['ASP', 'GLU', 'HIS']

superData = makeSuperDict([sims, []])
macroPkas = {'ASP': 3.65, 'GLU': 4.25, 'HIS': 6.53}

for sim in sims:
    if sim in ['4HFI_4', '6ZGD_4']:
        pH = 4.0
    else:
        pH = 7.0
    for resname in resnames:
        superData[sim].append(round(theoreticalProtonation(pH, macroPkas[resname]), 2))

ncols = 3
nrows = 1
data = resnames
outname = 'theoretical'
fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(1.8 * ncols, 4.5 * nrows), dpi=300)

row = 0
col = 0

for idx in range(0, len(data)):

    # Select subplt object.
    if nrows > 1:
        subplt = axs[row, col]
    else:
        subplt = axs[col]

    # The (recurring) bar plot block.
    width = 0.3
    x = np.arange(len([1]))

    mean4 = [superData['6ZGD_7'][idx]]
    mean3 = [superData['6ZGD_4'][idx]]
    mean2 = [superData['4HFI_7'][idx]]
    mean1 = [superData['4HFI_4'][idx]]

    #? ALL FOUR
    # subplt.bar(x - width * 1.5, mean4, width, color='#8856a7')
    # subplt.bar(x - width / 2.0, mean3, width, color='#9ebcda')
    # subplt.bar(x + width / 2.0, mean2, width, color='#8856a7', edgecolor='w', lw=1, hatch='//')
    # subplt.bar(x + width * 1.5, mean1, width, color='#9ebcda', edgecolor='w', lw=1, hatch='\\\\')

    #? JUST TWO
    subplt.bar(x - 0.5 * width, mean4, width, color='#8856a7')
    subplt.text(x - 0.5 * width - 0.14, mean4[0] + 0.01, mean4[0], size=15)

    subplt.bar(x + 0.5 * width, mean3, width, color='#9ebcda')
    subplt.text(x + 0.5 * width - 0.14, mean3[0] + 0.01, mean3[0], size=15)

    # Disable x-ticks and set x-label.
    subplt.set_xticks([])
    letter = [resnames[idx]]
    subplt.set_xlabel(['ASP', 'GLU', 'HIS'][idx])

    # Set y-lim.
    subplt.set_ylim([0, 1.1])

    # Remove boxes around plots.
    for loc in ['top', 'right', 'bottom', 'left']:
        subplt.spines[loc].set_visible(False)

    # If we're not in the first column, do not show the yticks.
    if col != 0:
        subplt.set_yticks([])
    else:
        subplt.set_ylabel('Protonation')
        subplt.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])

    # Increment row and column indices.
    col += 1
    if col == ncols:
        row += 1
        col = 0

#* REMOVE EMPTY BOXES IF LAST ROW WASN'T FILLED.
if (nrows > 1) and (len(data) % ncols != 0):

    left = ncols - len(data) % ncols
    # print(len(data), ncols, left)  # debug

    for ii in range(0, left):
        axs[nrows - 1, ncols - 1 - ii].axis('off')

#* LEGEND
subplt.bar(0, 0, 0, color='#8856a7', label='pH 7.0')
subplt.bar(0, 0, 0, color='#9ebcda', label='pH 4.0')
# subplt.bar(0, 0, 0, color='#8856a7', label='open, pH 7', edgecolor='w', lw=1, hatch='//')
# subplt.bar(0, 0, 0, color='#9ebcda', label='open, pH 4', edgecolor='w', lw=1, hatch='\\\\')
fig.legend(loc=[0.25, 0.65])

#* SAVE
# fig.legend()
fig.tight_layout()
fig.savefig(f'{outname}.png')
fig.savefig(f'{outname}.eps')
fig.clf()
