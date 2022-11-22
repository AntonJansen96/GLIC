import matplotlib.pyplot as plt
import numpy as np
import MDAnalysis as mda


def simple_dist(coord1, coord2):
    return np.sqrt((coord1[0] - coord2[0])**2 + (coord1[1] - coord2[1])**2 + (coord1[2] - coord2[2])**2)


dataset = '4HFI_4_1'  # not sure whether this variable is necessary
top     = '../../sims/4HFI_4/01/CA.pdb'
trajs   = ['../../sims/4HFI_4/01/MD_conv.xtc']
SAVE    = True

ecd_twist            = []
ecd_spread           = []
ecd_upper_spread     = []
beta_expansion       = []
m2_m1_dist           = []
m2_radius            = []
m1_kink              = []
m1_kink_alt          = []
nine_prime_dist      = []
nine_prime_pore      = []
minus_two_prime_dist = []
minus_two_prime_pore = []
c_loop               = []

for traj in trajs:

    ecd_twist1            = [[[traj, 0]], [[traj, 1]], [[traj, 2]], [[traj, 3]], [[traj, 4]]]
    ecd_spread1           = [[[traj, 0]], [[traj, 1]], [[traj, 2]], [[traj, 3]], [[traj, 4]]]
    ecd_upper_spread1     = [[[traj, 0]], [[traj, 1]], [[traj, 2]], [[traj, 3]], [[traj, 4]]]
    beta_expansion1       = [[[traj, 0]], [[traj, 1]], [[traj, 2]], [[traj, 3]], [[traj, 4]]]
    m2_m1_dist1           = [[[traj, 0]], [[traj, 1]], [[traj, 2]], [[traj, 3]], [[traj, 4]]]
    m2_radius1            = [[[traj, 0]], [[traj, 1]], [[traj, 2]], [[traj, 3]], [[traj, 4]]]
    m1_kink1              = [[[traj, 0]], [[traj, 1]], [[traj, 2]], [[traj, 3]], [[traj, 4]]]
    m1_kink_alt1          = [[[traj, 0]], [[traj, 1]], [[traj, 2]], [[traj, 3]], [[traj, 4]]]
    nine_prime_dist1      = [[[traj, 0]], [[traj, 1]], [[traj, 2]], [[traj, 3]], [[traj, 4]]]
    nine_prime_pore1      = [[[traj, 0]], [[traj, 1]], [[traj, 2]], [[traj, 3]], [[traj, 4]]]
    minus_two_prime_dist1 = [[[traj, 0]], [[traj, 1]], [[traj, 2]], [[traj, 3]], [[traj, 4]]]
    minus_two_prime_pore1 = [[[traj, 0]], [[traj, 1]], [[traj, 2]], [[traj, 3]], [[traj, 4]]]
    c_loop1               = [[[traj, 0]], [[traj, 1]], [[traj, 2]], [[traj, 3]], [[traj, 4]]]

    # Create universe
    u = mda.Universe(top, traj)
    nframes = u.trajectory.n_frames

    # Selections
    CA = u.select_atoms("chainID A B C D E and name CA")
    ECD_resids = "resid 0:193"
    TMD_resids = "resid 194:315"
    ECD_upper_resids = CA.select_atoms("resid 0:22") + CA.select_atoms("resid 43:78") + CA.select_atoms("resid 85:106") + CA.select_atoms("resid 131:149") + CA.select_atoms("resid 170:184")

    # Divide into chains
    previous = 0
    i = 0
    chains = [[], [], [], [], []]
    for ca in list(CA):
        current = ca.resid
        if current < previous:
            i += 1
        chains[i].append(ca)
        previous = current

    # Make each chain an AtomGroup. chains is not now a list of 5 AtomGroup
    # objects (one for each subunit), each containing 311 (CA) atoms.
    chains = [mda.AtomGroup(chains[i]) for i in range(0, len(chains))]

    #? SO FAR SO GOOD

    for j in range(0, nframes):

        ECD_com             = CA.select_atoms(ECD_resids).center_of_mass()
        TMD_com             = CA.select_atoms(TMD_resids).center_of_mass()
        ECD_upper_com       = ECD_upper_resids.center_of_mass()
        nine_prime_COM      = CA.select_atoms("resid 233").center_of_mass()
        minus_two_prime_COM = CA.select_atoms("resid 222").center_of_mass()

        for i in range(0, len(chains)):
            # Calculate ECD twist
            ECD_su_com = chains[i].select_atoms(ECD_resids).center_of_mass()
            TMD_su_com = chains[i].select_atoms(TMD_resids).center_of_mass()
            ECD_su_upper_com = (chains[i].select_atoms("resid 0:22") + chains[i].select_atoms("resid 43:78") + chains[i].select_atoms("resid 85:106") + chains[i].select_atoms("resid 131:149") + chains[i].select_atoms("resid 170:184")).center_of_mass()
            ecd_twist_coords = np.array([ECD_su_com, ECD_com, TMD_com, TMD_su_com])
            ecd_twist_universe = mda.Universe.empty(4, trajectory=True)
            ecd_twist_universe.atoms.positions = ecd_twist_coords
            ecd_twist1[i].append(mda.core.topologyobjects.Dihedral([0, 1, 2, 3], ecd_twist_universe).dihedral())

            # Calculate ECD spread
            ecd_spread1[i].append(simple_dist(ECD_com, ECD_su_com))

            # Calculate ECD upper spread
            ecd_upper_spread1[i].append(simple_dist(ECD_upper_com, ECD_su_upper_com))

            # Calculate beta expansion
            beta_com1 = chains[i].select_atoms("resid 30:34").center_of_mass()
            beta_com2 = chains[i].select_atoms("resid 190:194").center_of_mass()
            beta_expansion1[i].append(simple_dist(beta_com1, beta_com2))

            # Calculate M2-M1 distance
            m1_com = chains[i].select_atoms("resid 200:204").center_of_mass()
            if i == len(chains) - 1:
                m2_com = chains[0].select_atoms("resid 241:245").center_of_mass()
            else:
                m2_com = chains[i + 1].select_atoms("resid 241:245").center_of_mass()
            m2_m1_dist1[i].append(simple_dist(m1_com, m2_com))

            # Calculate M2 radius
            m2_com = chains[i].select_atoms("resid 241:245").center_of_mass()
            m2_TMD_com = CA.select_atoms("resid 241:245").center_of_mass()
            # M2 radius
            #m2_radius1[i].append(simple_dist(m2_com, m2_TMD_com))
            # M2 spread
            m2_radius1[i].append(simple_dist(m2_com, TMD_com))

            # Calculate M1 kink (apex at atom1)
            m1_com2 = chains[i].select_atoms("resid 202:206").center_of_mass()
            m1_com3 = chains[i].select_atoms("resid 205:216").center_of_mass()
            m1_kink_coords = np.array([m1_com, m1_com2, m1_com3])
            m1_kink_universe = mda.Universe.empty(3, trajectory=True)
            m1_kink_universe.atoms.positions = m1_kink_coords
            m1_kink1[i].append(mda.core.topologyobjects.Angle([0, 1, 2], m1_kink_universe).angle())

            # Alternative M1 kink calculation (top COM higher up)
            m1_com_alt = chains[i].select_atoms("resid 196:202").center_of_mass()
            m1_kink_alt_universe = mda.Universe.empty(3, trajectory=True)
            m1_kink_alt_universe.atoms.positions = np.array([m1_com_alt, m1_com, m1_com3])
            m1_kink_alt1[i].append(mda.core.topologyobjects.Angle([0, 1, 2], m1_kink_alt_universe).angle())

            # Calculate 9' distance
            # Note: center of mass is for CA (will be same as position of CA, but better format)
            nine_prime1 = chains[i].select_atoms("resid 233").center_of_mass()
            if i >= 3:
                nine_prime2 = chains[i - 3].select_atoms("resid 233").center_of_mass()
            else:
                nine_prime2 = chains[i + 2].select_atoms("resid 233").center_of_mass()
            nine_prime_dist1[i].append(simple_dist(nine_prime1, nine_prime2))

            # Alternative 9' calculation
            nine_prime_pore1[i].append(simple_dist(nine_prime1, nine_prime_COM))

            # Calculate -2' distance
            minus_two_prime1 = chains[i].select_atoms("resid 222").center_of_mass()
            if i >= 3:
                minus_two_prime2 = chains[i - 3].select_atoms("resid 222").center_of_mass()
            else:
                minus_two_prime2 = chains[i + 2].select_atoms("resid 222").center_of_mass()
            minus_two_prime_dist1[i].append(simple_dist(minus_two_prime1, minus_two_prime2))

            # Alternative -2' calculation
            minus_two_prime_pore1[i].append(simple_dist(minus_two_prime1, minus_two_prime_COM))

            # Calculate C-loop distance
            cloop_com = chains[i].select_atoms("resid 174:182").center_of_mass()
            if i == 0:
                turn_com = chains[(len(chains) - 1)].select_atoms("resid 145:150").center_of_mass()
            else:
                turn_com = chains[i - 1].select_atoms("resid 145:150").center_of_mass()
            c_loop1[i].append(simple_dist(cloop_com, turn_com))

        if j < nframes - 1:
            u.trajectory.next()

    for k in range(0, len(chains)):
        ecd_twist.append(ecd_twist1[k])
        ecd_spread.append(ecd_spread1[k])
        ecd_upper_spread.append(ecd_upper_spread1[k])
        beta_expansion.append(beta_expansion1[k])
        m2_m1_dist.append(m2_m1_dist1[k])
        m2_radius.append(m2_radius1[k])
        m1_kink.append(m1_kink1[k])
        m1_kink_alt.append(m1_kink_alt1[k])
        nine_prime_dist.append(nine_prime_dist1[k])
        nine_prime_pore.append(nine_prime_pore1[k])
        minus_two_prime_dist.append(minus_two_prime_dist1[k])
        minus_two_prime_pore.append(minus_two_prime_pore1[k])
        c_loop.append(c_loop1[k])

    # print(len(m2_radius))
    # print(len(m2_radius[0]))
    # print(len(m2_radius[1]))
    # print(m2_radius[0][1:])
    # print(m2_radius[0][0])
    # print(len(m2_radius[0][1:]))

    chains = ['A', 'B', 'C', 'D', 'E']
    for idx in range(0, len(chains)):

        t = range(0, len(m2_radius[idx][1:]))
        x = m2_radius[idx][1:]
        plt.plot(t, x, label=chains[idx])

    plt.legend()
    plt.xlabel('Time (ns)')
    plt.ylabel('Distance (A??)')
    plt.title('4HFI_4 1 chain A m2_radius')
    plt.tight_layout()
    plt.savefig('m2_radius_test.png')

# if SAVE:
#     ecd_twist = np.array(ecd_twist)
#     ecd_spread = np.array(ecd_spread)
#     ecd_upper_spread = np.array(ecd_upper_spread)
#     beta_expansion = np.array(beta_expansion)
#     m2_m1_dist = np.array(m2_m1_dist)
#     m2_radius = np.array(m2_radius)
#     m1_kink = np.array(m1_kink)
#     m1_kink_alt = np.array(m1_kink_alt)
#     nine_prime_dist = np.array(nine_prime_dist)
#     nine_prime_pore = np.array(nine_prime_pore)
#     minus_two_prime_dist = np.array(minus_two_prime_dist)
#     minus_two_prime_pore = np.array(minus_two_prime_pore)
#     c_loop = np.array(c_loop)

#     print(m2_radius)

#     '''
#     np.save("ecd_twist_" + dataset + ".npy", ecd_twist)
#     np.save("ecd_spread_" + dataset + ".npy", ecd_spread)
#     np.save("ecd_upper_spread_" + dataset + ".npy", ecd_upper_spread)
#     np.save("beta_expansion_" + dataset + ".npy", beta_expansion)
#     np.save("m2_m1_dist_" + dataset + ".npy", m2_m1_dist)
#     np.save("m2_radius_" + dataset + ".npy", m2_radius)
#     np.save("m1_kink_" + dataset + ".npy", m1_kink)
#     np.save("m1_kink_alt_" + dataset + ".npy", m1_kink_alt)
#     np.save("nine_prime_" + dataset + ".npy", nine_prime_dist)
#     np.save("nine_prime_pore_" + dataset + ".npy", nine_prime_pore)
#     np.save("minus_two_prime_" + dataset + ".npy", minus_two_prime_dist)
#     np.save("minus_two_prime_pore_" + dataset + ".npy", minus_two_prime_pore)
#     np.save("c_loop_" + dataset + ".npy", c_loop)
#     '''
#     np.save("m2_spread_" + dataset + ".npy", m2_radius)

    # print(list(m2_radius[1]))
    # print(len(list(m2_radius[1])))

    # plt.plot(range(0, len(list(m2_radius[1]))), list(m2_radius[1]))
    # plt.savefig('test.png')
