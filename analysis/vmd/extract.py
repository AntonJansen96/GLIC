import os
import MDAnalysis

for sim in ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']:

    os.chdir(sim)

    u   = MDAnalysis.Universe('{}.pdb'.format(sim))
    sel = u.select_atoms("(resname ASPT and name HD2) or (resname GLUT and name HE2)")

    with open('occupancies.txt', 'w+') as file:
        for ts in u.trajectory:
            for idx in list(sel.indices):
                file.write("{:.2f}\n".format(float(ts.data['occupancy'][idx])))

    os.chdir('..')
