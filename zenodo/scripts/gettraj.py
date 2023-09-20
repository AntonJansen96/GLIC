import os
import multiprocessing as mp

sims = ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']
reps = [1, 2, 3, 4]

# Get the .pdb and .xtc files.
for sim in sims:
    os.system(f"cp /home/anton/GIT/GLIC/sims/{sim}/01/CA.pdb ../production_traj/{sim}.pdb")
    for rep in reps:
        p1 = f"/home/anton/GIT/GLIC/sims/{sim}/{rep:02d}/MD_conv.xtc"
        p2 = f"../production_traj/{sim}_{rep}.xtc"
        os.system(f"gmx trjconv -f {p1} -o {p2} -dt 10000")

        os.mkdir(f"../production_traj/{sim}_{rep}_cph")

# Get the lambda coordinate files.
def task(sim, rep, coord):
    print(sim, rep, coord)

    source = open(f"/home/anton/GIT/GLIC/sims/{sim}/{rep:02d}/cphmd-coord-{coord}.xvg")
    target = open(f"../production_traj/{sim}_{rep}_cph/cphmd-coord-{coord}.xvg", 'w+')

    for line in source.readlines():
        if line[0] in ['#', '@']:
            target.write(line)
            continue

        num = float(line.split()[0])
        if num % 1000 == 0 and num <= 1e6:
            target.write(line)

    source.close()
    target.close()

items = []
for sim in sims:
    for rep in reps:
        for coord in range(1, 217):
            items.append((sim, rep, coord))

pool = mp.Pool(processes=mp.cpu_count())
pool.starmap(task, items, chunksize=1)
