# Software

* The GROMACS CpHMD beta can be found at www.gitlab.com/gromacs-constantph.

# analysis_results

* Analysis data and Python scripts for generating the (relevant) paper figures.

# prep_file

* Topology, force field, and GROMACS parameter files for preparing the simulations.
* Suggested to use the included updatelambdas.py to use lambda coordinates from previous EQ step as initial values for the next EQ step.

# production_traj

* 1 microsecond (10 ns/frame) trajectories for each of the four simulation systems and replicas.
* File names designations: 4HFI = open state structure, 6ZGD = closed state structure, second number is pH, third number is replica.
* The cph directories contain the 216 lambda coordinates.
  * 1 microsecond (1 ns/frame).
  * 5 subunits x 43 coordinates (1 for each ASPT/GLUT, 3 for each HSPT) + 1 BUF = 216.
