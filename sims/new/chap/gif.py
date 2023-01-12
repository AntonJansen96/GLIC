import matplotlib
import matplotlib.pyplot as plt
import multiprocessing as mp
import json
import os

from science.utility import gromacs
from science.utility import inputOptionHandler

matplotlib.rcParams.update({'font.size': 16})

options = ['Do nothing', 'Split trajectory', 'Jobscript for backend node', 'Jobscript for login node', 'Make images', 'Make gif']
opmode = inputOptionHandler('Select part of program to run', options)

# We can only do this for one replica at a time since we now have a time dependence.
# PART 1: Split the trajectories of new/4HFI_4/01 and new/4HFI_7/01 into 1ns pieces.

items = range(0, 301)  # First 300ns.

def task(frame):
    gromacs(f'trjconv -f ./../4HFI_4/01/MD_conv.xtc -o gif/4HFI_4_{frame}.xtc -b {1000 * frame} -e {1000 * frame}')
    gromacs(f'trjconv -f ./../4HFI_7/01/MD_conv.xtc -o gif/4HFI_7_{frame}.xtc -b {1000 * frame} -e {1000 * frame}')

if __name__ == "__main__":

    if opmode == 1:
        pool = mp.Pool(processes=mp.cpu_count())
        pool.map(task, items, chunksize=1)

# PART 2: Write a jobscript file for all.
################################################################################

if opmode == 2:
    with open('gif/jobs.txt', 'w+') as file:
        for frame in items:
            file.write(f'chap -s CA.pdb -f 4HFI_4_{frame}.xtc -n CHAP.ndx -out-filename 4HFI_4_{frame} -pf-sel-ipp 32 -sel-pathway 1 -pf-vdwr-fallback 0 -pf-max-free-dist 2 -pf-cutoff 2 -hydrophob-fallback 0\n')
            file.write(f'chap -s CA.pdb -f 4HFI_7_{frame}.xtc -n CHAP.ndx -out-filename 4HFI_7_{frame} -pf-sel-ipp 32 -sel-pathway 1 -pf-vdwr-fallback 0 -pf-max-free-dist 2 -pf-cutoff 2 -hydrophob-fallback 0\n\n')

# PART 2 : Alternative: run in parallel on login node.
################################################################################
# Use "ls job_*.sh | while read i; do bash $i & done" to start them.

if opmode == 3:
    frame = 0
    num = 1
    while True:
        file = open(f'gif/job_{num}.sh', 'w+')
        file.write('''module unload gromacs
    module load gromacs/2018.6
    module load lapack/3.8.0
    module load hwloc/1.11.13
    module load chap/0.9.1-gromacs2018.6\n\n''')
        count = 0

        while (count != 15):
            if frame > 300:
                break

            file.write(f'chap -s CA.pdb -f 4HFI_4_{frame}.xtc -n CHAP.ndx -out-filename 4HFI_4_{frame} -pf-sel-ipp 32 -sel-pathway 1 -pf-vdwr-fallback 0 -pf-max-free-dist 2 -pf-cutoff 2 -hydrophob-fallback 0\n')
            file.write(f'chap -s CA.pdb -f 4HFI_7_{frame}.xtc -n CHAP.ndx -out-filename 4HFI_7_{frame} -pf-sel-ipp 32 -sel-pathway 1 -pf-vdwr-fallback 0 -pf-max-free-dist 2 -pf-cutoff 2 -hydrophob-fallback 0\n\n')
            frame += 1
            count += 1

        num += 1
        if frame > 300:
            break

# PART 3: Make the images.
################################################################################

def swap(*line_list):
    for lines in line_list:
        try:
            iter(lines)
        except:  # noqa
            lines = [lines]
        for line in lines:
            xdata, ydata = line.get_xdata(), line.get_ydata()
            line.set_xdata(ydata)
            line.set_ydata(xdata)
            line.axes.autoscale_view()

def task(frame):
    plt.figure(figsize=(5, 7))

    data = json.load(open('4HFI_static/protein.json'))
    x = [i * 10 for i in data['pathwayProfile']['s']]
    y = [i * 10 for i in data['pathwayProfile']['radiusMean']]

    line1 = plt.plot(x, y, lw=1.5, label='4HFI crystal', color='black')
    swap(line1)

    for sim in ['4HFI_4', '4HFI_7']:

        data = json.load(open(f'gif/{sim}_{frame}.json'))

        x = [i * 10 for i in data['pathwayProfile']['s']]
        y = [i * 10 for i in data['pathwayProfile']['radiusMean']]

        if sim == '4HFI_4':
            line1 = plt.plot(x, y, lw=1.5, label=f'{sim} New', color='#1f77b4')
        elif sim == '4HFI_7':
            line1 = plt.plot(x, y, lw=1.5, label=f'{sim} New', color='#ff7f0e')
        swap(line1)

    plt.axis([0, 9, -20, 50])
    plt.grid()
    plt.xlabel('Radius (Å)')
    plt.ylabel('Z-coordinate (Å)')
    plt.legend(fontsize='small')
    plt.title('protein, frame = {:03d} ns'.format(frame))
    plt.tight_layout()
    plt.savefig('gif/{:03d}.png'.format(frame))
    plt.clf()
    plt.close()

if __name__ == "__main__":
    if opmode == 4:
        pool = mp.Pool(processes=mp.cpu_count())
        pool.map(task, items, chunksize=1)

# PART 5: Make the gif.
################################################################################

if opmode == 5:
    os.chdir('gif')
    os.system('ffmpeg -framerate 15 -pattern_type glob -i \'*.png\' -r 15 -vf scale=512:-1 chap1.gif')
    os.system('convert -delay 5 -loop 0 *.png chap2.gif')
