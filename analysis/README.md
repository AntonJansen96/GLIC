`analyze.py`
* Associated with the directory `lambdaplots`.
* Requires presence of `data.py` and `data.pickle` (optional).
* Makes the older histogram plots amongst other things. Still useful.

`contactData.py`
* Associated with the directory `contacts`.
* Performs systematic contact occupancy analysis.
* Creates the xxx.dat overview files as well as raw xxx.txt data.
* Multithreaded and requires a couple hours to run for all.

`contactPlot.py`
* Associated with the directory `contacts`.
* Requires prior running of `contactData.py` as it needs the `contacts/xxx` files.
* BUG: protonation columns for HSPT (in the plots) are not correct.

`data.pickle`
* Holds some data structures for speedup. Used by `analyze.py`.

`data.py`
* Experimental data on protonation states. Used by `analyze.py`.

`doPanels.py`
* Associated with the directory `panels`.
* Creates the older panels for analyzing residue contacts. Now mostly redundant.

`getNames.py`
* Helper file for identiying what kind of selections MDAnalysis makes.
* This was used to do some debugging for `contactData.py`.

`specialPlot.py`
* Associated with the directory `contacts`.
* Requires prior running of `timeData.py` as it needs the `time/xxx.xvg` files.
* Makes more advanced bar plots where we use different binning schemes.
* BUG: The entire plots for HSPT are incorrect.

`timeData.py`
* Associated with the directory `time`.
* Requires prior running of `contactData.py` as it needs the `contacts/xxx.dat`  files.
* Gathers the mindist data required for creating the in-time and special bar plots.
* Multithreaded and requires around 8 hours to run for all.

`timePlot.py`
* Associated with the directory `time`.
* Requires prior running of `timeData.py` as it needs the `time/xxx.xvg` files.
* Creates the 4x5 panels for better insight in proto-structure relationship in time.
* BUG: Protonation columns for HSPT are not correct.

chap
* Everything related to the CHAP analysis.

rmsd
* RMSD of the ECD (residues 18 to 200).
* Used early on in GLIC analysis process to check for simulation stability.

vmd
* Everything related to visual analysis and rendering.