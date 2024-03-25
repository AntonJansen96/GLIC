import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

fancy  = ["Closed, pH 7.0", "Closed, pH 4.0", "Open, pH 7.0", "Open, pH 4.0"]
sims   = ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']
reps   = [1, 2, 3, 4]
chains = ['A', 'B', 'C', 'D', 'E']

# THIS PART LOADS THE .CSV FILE ################################################

df = pd.read_csv('rmsd.csv')

# THIS PART CREATES SUPPLEMENTAL FIGURE 1 ######################################

matplotlib.rcParams.update({'font.size': 22})

nrows = 4
ncols = 4

fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(17, 11), dpi=200)

for row in range(0, nrows):
    for col in range(0, ncols):

        subplt = axs[row, col]

        for chain in chains:
            t = range(0, 1001)
            x = df[f"{sims[col]}_{reps[row]}_{chain}"]
            # x = loadCol(f"{sims[col]}_{reps[row]}_{chain}.txt", col=2)
            subplt.plot(t, x, linewidth=0.5)

        # Set x-lim and y-lim.
        subplt.set_ylim([0.5, 4.5])
        subplt.set_xlim([0, 1000])

        subplt.text(15, 4.0, str(row+1))

        # If we're not in the last row, do not show the xticks.
        if row != nrows - 1:
            subplt.set_xticks([])
        else:
            subplt.set_xlabel("Time (ns)")
            subplt.set_xticks([0, 250, 500, 750, 1000])

        # If we're not in the first column, do not show the yticks.
        if col != 0:
            subplt.set_yticks([])
        else:
            subplt.set_ylabel(r'ECD RMSD ($\AA$)', size=20)
            subplt.set_yticks([1, 2, 3, 4])

        # Add title to top row.
        if row == 0:
            subplt.set_title(fancy[col])

fig.tight_layout(pad=0.3)
fig.savefig('supplement1.png')
fig.clf()
