import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from functions import makeSuperDict

sims     = ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']
residues = [13, 14, 26, 31, 32, 35, 49, 55, 67, 69, 75, 82, 86, 88, 91, 97, 104, 115, 122, 127, 136, 145, 147, 153, 154, 161, 163, 177, 178, 181, 185, 222, 235, 243, 272, 277, 282]

# THIS PART LOADS THE .CSV FILE ################################################

protoMean = makeSuperDict([sims, residues, []])

df = pd.read_csv('protonation.csv')
for sim in sims:
    for res in residues:
            protoMean[sim][res] = list(df[f'{sim}_{res}'])

# THIS PART CREATES FIGURE 2 ###################################################

matplotlib.rcParams.update({'font.size': 24})

def stderr(array: list) -> float:
    return np.std(array) / np.sqrt(len(array))

data = ['D13', 'E14', 'E67', 'E69', 'E104', 'H127', 'E147', 'D154', 'D161', 'D178', 'E222']

protoErr = makeSuperDict([sims, residues, 0])

for sim in sims:
    for residue in residues:
        protoErr[sim][residue] = stderr(protoMean[sim][residue])
        protoMean[sim][residue] = np.mean(protoMean[sim][residue])

ncols = 6
nrows = 2

fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(1.8 * ncols, 4.5 * nrows), dpi=300)

row = 0
col = 0

for idx in range(0, len(data)):

    subplt = axs[row, col]

    width = 0.3
    x = np.arange(len([1]))

    if data[idx] in ['ASP', 'GLU', 'HIS']:
        mean1 = [{'ASP': 0.03, 'GLU': 0.06, 'HIS': 0.38}[data[idx]]]
        mean2 = [{'ASP': 0.41, 'GLU': 0.56, 'HIS': 0.93}[data[idx]]]

        subplt.bar(x - 0.5 * width, mean1, width, color='#8856a7')
        subplt.text(x - 0.5 * width - 0.14, mean1[0] + 0.01, mean1[0], size=19)

        subplt.bar(x + 0.5 * width, mean2, width, color='#9ebcda')
        subplt.text(x + 0.5 * width - 0.14, mean2[0] + 0.01, mean2[0], size=19)

    else:
        mean4 = [protoMean['6ZGD_7'][float(data[idx][1:])]]
        serr4 = [protoErr['6ZGD_7'][float(data[idx][1:])]]
        subplt.bar(     x - width * 1.5, mean4, width, color='#8856a7')
        subplt.errorbar(x - width * 1.5, mean4, serr4, color='#8856a7', fmt='none', capsize=7, linewidth=3, markeredgewidth=2)

        mean3 = [protoMean['6ZGD_4'][float(data[idx][1:])]]
        serr3 = [protoErr['6ZGD_4'][float(data[idx][1:])]]
        subplt.bar(     x - width / 2.0, mean3, width, color='#9ebcda')
        subplt.errorbar(x - width / 2.0, mean3, serr3, color='#9ebcda', fmt='none', capsize=7, linewidth=3, markeredgewidth=2)

        mean2 = [protoMean['4HFI_7'][float(data[idx][1:])]]
        serr2 = [protoErr['4HFI_7'][float(data[idx][1:])]]
        subplt.bar(     x + width / 2.0, mean2, width, color='#8856a7', edgecolor='w', lw=0, hatch='//', zorder=13)
        subplt.errorbar(x + width / 2.0, mean2, serr2, color='#8856a7', fmt='none', capsize=7, linewidth=3, markeredgewidth=2, zorder=10)

        mean1 = [protoMean['4HFI_4'][float(data[idx][1:])]]
        serr1 = [protoErr['4HFI_4'][float(data[idx][1:])]]
        subplt.bar(     x + width * 1.5, mean1, width, color='#9ebcda', edgecolor='w', lw=0, hatch='\\\\', zorder=12)
        subplt.errorbar(x + width * 1.5, mean1, serr1, color='#9ebcda', fmt='none', capsize=7, linewidth=3, markeredgewidth=2, zorder=11)

    subplt.set_xticks([])           # Disable x-ticks and set x-label.
    subplt.set_xlabel(data[idx])
    subplt.set_ylim([0, 1.0])       # Set y-lim.

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

# Remove empty boxes if the last row wasn't filled.
if (nrows > 1) and (len(data) % ncols != 0):
    left = ncols - len(data) % ncols

    for ii in range(0, left):
        axs[nrows - 1, ncols - 1 - ii].axis('off')

fig.tight_layout(pad=0.8)
fig.savefig('supplement2.png')
fig.clf()
