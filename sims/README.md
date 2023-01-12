* Make dirs, copy POPC.itp, run.sh, equilibrate.sh, jobscript.sh

For each:

* phbuilder gentopol -f 4HFI_charmmgui.pdb -auto xxx
* -ter
* neutral terminals (option 2 for both NTER and CTER)

* change name TIP3 to SOL in phprocessed.pdb
* change name TIP3 in topol.top (only at the bottom)

* phbuilder neutralize -f phprocessed.pdb -conc 0.15
* phbuilder genparams -f phneutral.pdb -ph xxx
* membrane protein = yes

* Change in NVT.mdp the steps and dt with 10
* run ./equilibrate.sh

* Change simulation time in MD.mdp to 20ns
* Update jobscript name

--------------------------------------------------------------------------------

* Get which residues' Edwp needs to be updated, and by how much
* Copy Edwp sims to a different folder
* phbuilder genparams -f phneutral.pdb -ph xxx -inter
* Change simulation time to 200ns, change nstout for .xtc
* cleangromacsbackups; rm slurm* MD.xtc

--------------------------------------------------------------------------------

* mkdir bck_1; cp lambda_*.dat ./bck_1/
* for i in {1..216}; do sed -i '$d' lambda_$i.dat; done
* (Edit the jobscript to do a restart with appending)
* sbatch jobscript

--------------------------------------------------------------------------------

For each sim
* rm MD.pdb (this is an old file from Edwpruns)
* gmx convert-tpr -s MD.tpr -o MD.tpr -extend 800000
* Update jobscript

* runBackup
* on tcb copy GLIC from runBackup to scratch on /cephmfs
* sbatch the jobscripts

--------------------------------------------------------------------------------

ANALYSIS (at least for update Feb 20th)

* Remove last line
* cleangromacsbackups
* rm -rf bin include lib share bck_1 bck_2 slurm*
* unset DISPLAY
* ./doanalysis.py

--------------------------------------------------------------------------------

Opnieuw opstarten
* gmx genrestr om restraint files te maken
* ifdefs toevoegen aan verschillen topol_xx.itp

Test
* Remove lipid restraint from NVT.mdp           yes
* Add -r xxx.pdb for all grompp lines           yes
* Divide time step by ten in NVT.mdp            
* Remove excess lines in BACK.itp and CA.itp    yes
* 