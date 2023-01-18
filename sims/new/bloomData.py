#!/bin/python3

import MDAnalysis as mda
import numpy as np
import multiprocessing as mp
import os


def doBlooming(sim, rep):
    """Run Cathrine's blooming analysis for a simulation."""

    def dist(c1, c2) -> float:
        """Euclidian distance between two points in 3D space."""
        x = c1[0] - c2[0]
        y = c1[1] - c2[1]
        z = c1[2] - c2[2]
        return (x * x + y * y + z * z)**0.5

    # Prepare directory stuff
    if not os.path.exists('bloom'):
        os.mkdir('bloom')

    # Get the correct pdb and xtc file.
    pdb = "{}/{:02d}/CA.pdb".format(sim, rep)
    xtc = "{}/{:02d}/MD_conv.xtc".format(sim, rep)
    u   = mda.Universe(pdb, xtc)

    # This is to get values for static structures. Normally commented out.
    # pdb = "../sims/source/4HFI_clean_crystal.pdb"
    # pdb = "../sims/source/6ZGD_clean_cryoEM.pdb"
    # u  = mda.Universe(pdb)

    # Here we define the analysis metrics and data structure for our
    # resulting data. Each metric is a list of 5 lists. Each list holds
    # N-frames of values for a specific chain.
    metrics   = ['ecd_twist', 'ecd_spread', 'ecd_upper_spread', 'beta_expansion', 'm2_m1_dist', 'm2_radius', 'm1_kink', 'm1_kink_alt', 'nine_prime_dist', 'nine_prime_pore', 'minus_two_prime_dist', 'minus_two_prime_pore', 'c_loop']
    superData = {}
    for key in metrics:
        superData[key] = [[], [], [], [], []]

    # Some selections
    CA = u.select_atoms("chainID A B C D E and name CA")
    ECD_resids = "resid 0:193"
    TMD_resids = "resid 194:315"
    ECD_upper_resids = CA.select_atoms("resid 0:22") + CA.select_atoms("resid 43:78") + CA.select_atoms("resid 85:106") + CA.select_atoms("resid 131:149") + CA.select_atoms("resid 170:184")

    # chains contains five mda arrays. Each contains Calpha atoms for one chain.
    chains = []
    for chain in ['A', 'B', 'C', 'D', 'E']:
        chains.append(u.select_atoms('chainID {} and name CA'.format(chain)))

    # Normally commented out: this is to get only the last 500ns.
    frameCount = 0
    for _ in u.trajectory:

        frameCount += 1
        if frameCount > 300:
            continue

        #? THIS PART I COPIED FROM CATHRINES SCRIPT AND I DON'T FULLY
        #? UNDERSTAND EVERY ANALYSIS METRIC.

        ECD_com             = CA.select_atoms(ECD_resids).center_of_mass()
        TMD_com             = CA.select_atoms(TMD_resids).center_of_mass()
        ECD_upper_com       = ECD_upper_resids.center_of_mass()
        nine_prime_COM      = CA.select_atoms("resid 233").center_of_mass()
        minus_two_prime_COM = CA.select_atoms("resid 222").center_of_mass()

        ca_com = u.select_atoms('resid 233 and name CA').center_of_mass()  #! Used for corrected nine_prime_dist computation.

        for i in range(0, len(chains)):
            # Calculate ECD twist
            ECD_su_com = chains[i].select_atoms(ECD_resids).center_of_mass()
            TMD_su_com = chains[i].select_atoms(TMD_resids).center_of_mass()
            ECD_su_upper_com = (chains[i].select_atoms("resid 0:22") + chains[i].select_atoms("resid 43:78") + chains[i].select_atoms("resid 85:106") + chains[i].select_atoms("resid 131:149") + chains[i].select_atoms("resid 170:184")).center_of_mass()
            ecd_twist_coords = np.array([ECD_su_com, ECD_com, TMD_com, TMD_su_com])
            ecd_twist_universe = mda.Universe.empty(4, trajectory=True)
            ecd_twist_universe.atoms.positions = ecd_twist_coords
            superData['ecd_twist'][i].append(mda.core.topologyobjects.Dihedral([0, 1, 2, 3], ecd_twist_universe).dihedral())

            # Calculate ECD spread
            superData['ecd_spread'][i].append(dist(ECD_com, ECD_su_com))

            # Calculate ECD upper spread
            superData['ecd_upper_spread'][i].append(dist(ECD_upper_com, ECD_su_upper_com))

            # Calculate beta expansion
            beta_com1 = chains[i].select_atoms("resid 30:34").center_of_mass()
            beta_com2 = chains[i].select_atoms("resid 190:194").center_of_mass()
            superData['beta_expansion'][i].append(dist(beta_com1, beta_com2))

            # Calculate M2-M1 distance
            m1_com = chains[i].select_atoms("resid 200:204").center_of_mass()
            if i == len(chains) - 1:
                m2_com = chains[0].select_atoms("resid 241:245").center_of_mass()
            else:
                m2_com = chains[i + 1].select_atoms("resid 241:245").center_of_mass()
            superData['m2_m1_dist'][i].append(dist(m1_com, m2_com))

            # Calculate M2 radius
            m2_com = chains[i].select_atoms("resid 241:245").center_of_mass()
            # m2_TMD_com = CA.select_atoms("resid 241:245").center_of_mass()
            # M2 radius
            #m2_radius1[i].append(dist(m2_com, m2_TMD_com))
            # M2 spread
            superData['m2_radius'][i].append(dist(m2_com, TMD_com))

            # Calculate M1 kink (apex at atom1)
            m1_com2 = chains[i].select_atoms("resid 202:206").center_of_mass()
            m1_com3 = chains[i].select_atoms("resid 205:216").center_of_mass()
            m1_kink_coords = np.array([m1_com, m1_com2, m1_com3])
            m1_kink_universe = mda.Universe.empty(3, trajectory=True)
            m1_kink_universe.atoms.positions = m1_kink_coords
            superData['m1_kink'][i].append(mda.core.topologyobjects.Angle([0, 1, 2], m1_kink_universe).angle())

            # Alternative M1 kink calculation (top COM higher up)
            m1_com_alt = chains[i].select_atoms("resid 196:202").center_of_mass()
            m1_kink_alt_universe = mda.Universe.empty(3, trajectory=True)
            m1_kink_alt_universe.atoms.positions = np.array([m1_com_alt, m1_com, m1_com3])
            superData['m1_kink_alt'][i].append(mda.core.topologyobjects.Angle([0, 1, 2], m1_kink_alt_universe).angle())

            # Calculate 9' distance
            # Note: center of mass is for CA (will be same as position of CA, but better format)
            nine_prime1 = chains[i].select_atoms("resid 233").center_of_mass()
            if i >= 3:
                nine_prime2 = chains[i - 3].select_atoms("resid 233").center_of_mass()
            else:
                nine_prime2 = chains[i + 2].select_atoms("resid 233").center_of_mass()
            # superData['nine_prime_dist'][i].append(dist(nine_prime1, nine_prime2))  #! Comments this out.

            # Alternative 9' calculation
            superData['nine_prime_pore'][i].append(dist(nine_prime1, nine_prime_COM))

            #! Anton: This comes from Cathrine's new scripts and fixes my problem
            #! with the nine_prime_dist metric (4HFI now gives 5.07 A, which is correct).
            letters = ['A', 'B', 'C', 'D', 'E']
            # Calculate -2' smallest distance
            min_dist = 10000000
            resid1 = u.select_atoms(f'protein and chainID {letters[i]} and resid 233')
            for pos in resid1.positions:
                distance = dist(pos, ca_com)
                if distance < min_dist:
                    min_dist = distance
            superData['nine_prime_dist'][i].append(min_dist)
            #! End of modification.

            # Calculate -2' distance
            minus_two_prime1 = chains[i].select_atoms("resid 222").center_of_mass()
            if i >= 3:
                minus_two_prime2 = chains[i - 3].select_atoms("resid 222").center_of_mass()
            else:
                minus_two_prime2 = chains[i + 2].select_atoms("resid 222").center_of_mass()
            superData['minus_two_prime_dist'][i].append(dist(minus_two_prime1, minus_two_prime2))

            # Alternative -2' calculation
            superData['minus_two_prime_pore'][i].append(dist(minus_two_prime1, minus_two_prime_COM))

            # Calculate C-loop distance
            cloop_com = chains[i].select_atoms("resid 174:182").center_of_mass()
            if i == 0:
                turn_com = chains[(len(chains) - 1)].select_atoms("resid 145:150").center_of_mass()
            else:
                turn_com = chains[i - 1].select_atoms("resid 145:150").center_of_mass()
            superData['c_loop'][i].append(dist(cloop_com, turn_com))

        #? END OF COPIED PART.

    # Write everything to file(s).
    chains = ['A', 'B', 'C', 'D', 'E']
    for chain in chains:
        with open('bloom/{}_{}_{}.txt'.format(sim, rep, chain), 'w') as file:
            # Write header.
            for metric in metrics:
                file.write('{} '.format(metric))
            file.write('\n')
            # Write data.
            for idx in range(0, len(superData['c_loop'][0])):
                for metric in metrics:
                    file.write('{:.4f} '.format(superData[metric][chains.index(chain)][idx]))
                file.write('\n')


if __name__ == "__main__":
    # GATHER ITERABLES
    items = []
    for sim in ['4HFI_4', '4HFI_7', '6ZGD_4', '6ZGD_7']:
        for rep in [1, 2, 3, 4]:
            items.append((sim, rep))

    # RUN MULTITHREADED
    pool = mp.Pool(processes=mp.cpu_count())
    pool.starmap(doBlooming, items, chunksize=1)

# This is to get values for static structures. Normally commented out.
# doBlooming('6ZGD', 0)
# doBlooming('4HFI', 0)
