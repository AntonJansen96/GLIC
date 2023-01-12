#!/bin/bash
#SBATCH --time=1-23:59:59
#SBATCH --nodes=1
#SBATCH -p lindahl3,lindahl4
#SBATCH --job-name=EQ_4HFI_4
#SBATCH --mail-user=anton.jansen@scilifelab.se
#SBATCH --mail-type=ALL
#SBATCH -G 1

# LOAD MODULES

module load cmake/latest
module load gcc/7.4
module load cuda/10.2

simdir="$PWD"

# COMPILE CUSTOM GROMACS ON SCRATCH

builddir="/scratch/$USER/$SLURM_JOBID/build"
mkdir -p "$builddir"
cd "$builddir"
CC=gcc-7 CXX=g++-7 cmake ~/gromacs-constantph -DGMX_USE_RDTSCP=ON -DCMAKE_INSTALL_PREFIX=${PWD}/.. -DGMX_BUILD_OWN_FFTW=ON -DGMX_GPU=CUDA
make -j 12
make install -j 12
cd ..
source ${PWD}/bin/GMXRC
gmxbasepath="${PWD}"

# RUN EQUILIBRATION BY CALLING smartEQ_cluster.py

cd "$simdir"
python3 smartEQ_cluster.py ${SLURM_JOB_CPUS_PER_NODE} ${gmxbasepath}
