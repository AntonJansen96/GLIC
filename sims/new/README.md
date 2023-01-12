General
* This directory was created on Tuesday December 20th 2022.
* Contains the 'new' equilibration and production runs for 4HFI_4, 4HFI_7, 6ZGD_4, 6ZGD_7.

Directory structure
* `ca_new.pdb`: a structure from Yuxuan, used to align and add the cavity waters.
* EQ contains the equilibration runs. This is only one per 4 replicas.
* 4HFI_4, 4HFI_7, 6ZGD_4, 6ZGD_7 contain the production runs (4 replicas each).

First, the following files were coppied from `sims/SIM/01` to `sims/new/EQ/SIM`:
* `EM.mdp`
* `NVT.mdp`
* `HEAVY.mdp`
* `BACK.mdp`
* `CA.mdp`
* `MD.mdp`
* `POPC.itp`
* `HEAVY.itp`
* `topol_Protein_chain_A.itp`
* `topol_Protein_chain_B.itp`
* `topol_Protein_chain_C.itp`
* `topol_Protein_chain_D.itp`
* `topol_Protein_chain_E.itp`
* `topol.top`
* `phneutral.pdb`
* `index.ndx`
* `phrecord.dat`
* `charmm36-mar2019-m6.ff`

Improvements with respect to the older EQ and production simulations:
* There was a bug in the position restraints of `BACK.itp` and `CA.itp`, namely that the atoms belonging
to the titratable residues (ASPT, GLUT, HSPT) were not in them. This can be problematic for EQ and has now been fixed.
* `MD.mdp` is now 1000ns instead of 2000ns, and `nstxout-compressed = 50000` (100ps) instead of `nstxout-compressed = 25000` (50ps) previously. This is to save drive space.
* In `EM.mdp`, `NVT.mdp`, `HEAVY.mdp`, `BACK.mdp`, `CA.mdp`, `nstxout-compressed = 50000` instead of `0`. </br>This allows for checking whether the restraints are actually working.
* The EQ has been made much more sound by using the lambda coordinate of the last frame as the starting lambda for the next EQ step. This is done by a function defined in `smartEQ.py`.
* A results of this is that we set `lambda-dynamics-calibration = no` for `NVT.mdp` and all subsequent steps.
* Another result of this is that the system is now much more stable and we were able to set `dt = 0.002` instead of `dt = 0.0002` for the `HEAVY.mdp` step.
* For some other EQ details, see EQ tab in GLIC mastersheet.

Furthermore, specifically for 4HFI_4 and 4HFI_7:
* Waters were added (in `phneutral.pdb`) in the cavity close to H235. This was done by aligning `ca_new.pdb` (provided by Yuxuan) and `phneutral.pdb` in VMD. In addition to `phneutral.pdb`, `topol.top` and `index.ndx` were also updated, as to account for the additional molecules in the system. The reason for adding the water molecules in the cavity is that it might prevent the pore collapse observed in earlier 4HFI simulations in CHARMM36.
