4HFI_raw_crystal.pdb    Is straight from RSCB.org
6ZGD_raw_cryoEM.pdb     Is straight from /mnt/cephfs/projects/2019050300_GLIC_gating_muts_MD/yzhuang/glic_cryo_paper_new/initial_models (renamed from ph7.pdb to ...)

4HFI_clean_crystal.pdb  Removed all the crap using addChainIdentifiers
6ZGD_clean_cryoEM.pdb   Removed all the crap using addChainIdentifiers

charmmgui parameters:

* Download source = OPM for both (with respective names)
* Orientation = Use PPM server
* Length of X and Y = 140, POPC = 1 : 1
* Deselect include ions
* Force Field Options = CHARMM36
* Input Generation Options = tick GROMACS
* Temperature = 300K

Resulting in:

* 4HFI_charmm-gui.tgz --> charmm-gui-4HFI
* 6ZGD_charmm-gui.tgz --> charmm-gui-6ZGD

Extract from these archives workable .pdb input files:

* gmx editconf -f step5_input.gro -o temp.pdb
* Get CRYST1 from temp.pdb and paste in step5_input.pdb
* Change TITLE
* addChainIdentifiers
* Change names of water atoms

Resulting in:

4HFI_charmmgui.pdb
6ZGD_charmmgui.pdb
