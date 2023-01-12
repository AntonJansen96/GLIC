#!/bin/bash

# Gromacs version to use:
source /usr/local/gromacs_constantph/bin/GMXRC

# Equilibration stages involved: EM, NVT, HEAVY, BACKBONE, CA, RIGID (optional)
# In rigid equilibration step, well-resolved residues were restrained.

# EM (100ps)
# gmx grompp -f EM.mdp -c phneutral.pdb -p topol.top -n index.ndx -o EM.tpr -maxwarn 1
# gmx mdrun -v -deffnm EM -c EM.pdb

# NVT COUPLING (100ps)
# gmx grompp -f NVT.mdp -c EM.pdb -p topol.top -n index.ndx -o NVT.tpr -r EM.pdb
# gmx mdrun -v -deffnm NVT -c NVT.pdb

# RESTRAIN HEAVY ATOMS (10ns)
gmx grompp -f HEAVY.mdp -c NVT.pdb -p topol.top -n index.ndx -o HEAVY.tpr -r NVT.pdb
gmx mdrun -v -deffnm HEAVY -c HEAVY.pdb

# RESTRAIN BACKBONE (10ns)
gmx grompp -f BACK.mdp -c HEAVY.pdb -p topol.top -n index.ndx -o BACK.tpr -r HEAVY.pdb
gmx mdrun -v -deffnm BACK -c BACK.pdb

# # RESTRAIN CA (10ns)
gmx grompp -f CA.mdp -c BACK.pdb -p topol.top -n index.ndx -o CA.tpr -r BACK.pdb
gmx mdrun -v -deffnm CA -c CA.pdb
