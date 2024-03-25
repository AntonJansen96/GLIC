import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from science.parsing import loadCol

fancy  = ["Closed, pH 7.0", "Closed, pH 4.0", "Open, pH 7.0", "Open, pH 4.0"]
sims   = ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']
reps   = [1, 2, 3, 4]
chains = ['A', 'B', 'C', 'D', 'E']

#! THIS PART CREATES THE .CSV FILE #############################################

df = pd.DataFrame()

for sim in sims:
    for rep in reps:
        for chain in chains:

            rmsd = loadCol(f"/home/anton/GIT/GLIC/analysis/revcom/{sim}_{rep}_{chain}.txt", col=2)

            if len(rmsd) == 999:                    # This is to take care of
                rmsd += [rmsd[998]] + [rmsd[998]]   # the fact that some sims
            elif len(rmsd) == 1000:                 # finished at 998 or 999ns
                rmsd += [rmsd[999]]                 # instead of 1000ns.

            df[f"{sim}_{rep}_{chain}"] = rmsd

df.to_csv('rmsd.csv')

#! THIS PART LOADS THE .CSV FILE ###############################################

df = pd.read_csv('rmsd.csv')

#! THIS PART CREATES SUPPLEMENTAL FIGURE 1 #####################################

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
        subplt.set_ylim([0, 7])
        subplt.set_xlim([0, 1000])

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
            subplt.set_ylabel(r'RMSD ($\AA$)')
            subplt.set_yticks([0, 2, 4, 6])

        # Add title to top row.
        if row == 0:
            subplt.set_title(fancy[col])

fig.tight_layout(pad=0.3)
fig.savefig('supplement1.png')
fig.clf()
