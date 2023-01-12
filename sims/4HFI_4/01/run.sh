#!/bin/bash

# Gromacs version to use:
source /usr/local/gromacs_constantph/bin/GMXRC

# gmx grompp -f MD.mdp -c CA.pdb -p topol.top -n index.ndx -o MD.tpr
gmx mdrun -v -deffnm MD -c MD.pdb -x MD.xtc -ntmpi 2 -npme 0 -
