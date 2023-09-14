#!/usr/bin/env python3

import multiprocessing as mp

def task(dir, coord):
    print(dir, coord)

    source = open(f"/home/anton/GIT/GLIC/sims/{dir}/01/cphmd-coord-{coord}.xvg", 'r')
    target = open(f"/home/anton/GIT/GLIC/zenodo/production_traj/{dir}_cphmd/cphmd-coord-{coord}.xvg", 'w+')

    for line in source.readlines():
        if line[0] in ['#', '@']:
            target.write(line)
            continue

        num = float(line.split()[0])
        if num % 1000 == 0 and num <= 1e6:
            target.write(line)

items = []
for dir in ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']:
    for coord in range(1, 217):
        items.append((dir, coord))

# Run multithreaded
pool = mp.Pool(processes=mp.cpu_count())
pool.starmap(task, items, chunksize=1)
