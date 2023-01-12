import os
import subprocess
import sys

cores = int(sys.argv[1])
basepath = sys.argv[2]


def gromacs(command, stdin=[], terminal=True, logFile='gromacs.log'):
    """Python function for handeling calls to GROMACS.

    Args:
        command (string): GROMACS command, e.g. 'make_ndx -f protein.pdb'.
        stdin (list, optional): list of input arguments. Defaults to [].
        terminal (bool, optional): print output to terminal. Defaults to True.
        basepath (str, optional): base path for GROMACS version to be used. Defaults to 'usr/local/gromacs_constantph'.
        logFile (str, optional): log file to use when terminal=False. Defaults to 'gromacs.log'.

    Returns:
        int: return code (0 if things were successful).
    """

    # If we don't pass envvars to subprocess (which happens by default) this works.
    path_to_gmx = os.path.normpath(basepath + '/' + 'bin/gmx')
    command = "{} {}".format(path_to_gmx, command)

    # Use the << EOF method to handle user input. Prepare string.
    if stdin:
        xstr = ' << EOF\n'
        for val in stdin:
            xstr += '{}\n'.format(val)
        command += xstr + 'EOF'

    if terminal:
        process = subprocess.run(command, shell=True, env={})
    else:
        with open(logFile, 'a+') as file:
            process = subprocess.run(command, shell=True, stdout=file, stderr=file, env={})

    if process.returncode != 0:
        if terminal:
            print("Failed to run \"{}\" (exitcode {}).".format(command, process.returncode))

        print("Failed to run \"{}\" (exitcode {}). Check your logfile ({})".format(command, process.returncode, logFile))

    return process.returncode


def extractAndUpdate(base, new):
    """Extract the last lambda-coordinate from a previous simulation (base) and
    use those as starting lambdas for the new (new) simulation.
    Args:
        base (str): old simulation (e.g. NVT).
        new (str): new simulation (e.g. NPT).
    """

    def lastLambdaVal(fname):
        "Retrieves the last lambda coordinate from a lambda file."

        with open(fname) as file:
            for line in file:
                pass
            return float(line.split()[1])

    # Directory stuff
    if not os.path.exists(base):
        os.mkdir(base)

    # Extract and move the lambda-coordinate files.
    gromacs('cphmd -s {}.tpr -e {}.edr -numplot 1'.format(base, base))
    os.system('mv cphmd-coord-*.xvg {}'.format(base))

    values = []
    for num in range(1, len(os.listdir(base)) + 1):
        val = lastLambdaVal('{}/cphmd-coord-{}.xvg'.format(base, num))
        values.append(val)
        print(val)  # debug

    acidics = ['ASPT', 'GLUT', 'GLUT', 'ASPT', 'ASPT', 'GLUT', 'ASPT', 'ASPT', 'GLUT', 'GLUT', 'GLUT', 'GLUT', 'ASPT', 'ASPT', 'ASPT', 'ASPT', 'GLUT', 'ASPT', 'ASPT', 'HSPT', 'ASPT', 'ASPT', 'GLUT', 'ASPT', 'ASPT', 'ASPT', 'GLUT', 'GLUT', 'ASPT', 'GLUT', 'ASPT', 'GLUT', 'HSPT', 'GLUT', 'GLUT', 'HSPT', 'GLUT', 'ASPT', 'GLUT', 'GLUT', 'ASPT', 'ASPT', 'GLUT', 'ASPT', 'ASPT', 'GLUT', 'GLUT', 'GLUT', 'GLUT', 'ASPT', 'ASPT', 'ASPT', 'ASPT', 'GLUT', 'ASPT', 'ASPT', 'HSPT', 'ASPT', 'ASPT', 'GLUT', 'ASPT', 'ASPT', 'ASPT', 'GLUT', 'GLUT', 'ASPT', 'GLUT', 'ASPT', 'GLUT', 'HSPT', 'GLUT', 'GLUT', 'HSPT', 'GLUT', 'ASPT', 'GLUT', 'GLUT', 'ASPT', 'ASPT', 'GLUT', 'ASPT', 'ASPT', 'GLUT', 'GLUT', 'GLUT', 'GLUT', 'ASPT', 'ASPT', 'ASPT', 'ASPT', 'GLUT', 'ASPT', 'ASPT', 'HSPT', 'ASPT', 'ASPT', 'GLUT', 'ASPT', 'ASPT', 'ASPT', 'GLUT', 'GLUT', 'ASPT', 'GLUT', 'ASPT', 'GLUT', 'HSPT', 'GLUT', 'GLUT', 'HSPT', 'GLUT', 'ASPT', 'GLUT', 'GLUT', 'ASPT', 'ASPT', 'GLUT', 'ASPT', 'ASPT', 'GLUT', 'GLUT', 'GLUT', 'GLUT', 'ASPT', 'ASPT', 'ASPT', 'ASPT', 'GLUT', 'ASPT', 'ASPT', 'HSPT', 'ASPT', 'ASPT', 'GLUT', 'ASPT', 'ASPT', 'ASPT', 'GLUT', 'GLUT', 'ASPT', 'GLUT', 'ASPT', 'GLUT', 'HSPT', 'GLUT', 'GLUT', 'HSPT', 'GLUT', 'ASPT', 'GLUT', 'GLUT', 'ASPT', 'ASPT', 'GLUT', 'ASPT', 'ASPT', 'GLUT', 'GLUT', 'GLUT', 'GLUT', 'ASPT', 'ASPT', 'ASPT', 'ASPT', 'GLUT', 'ASPT', 'ASPT', 'HSPT', 'ASPT', 'ASPT', 'GLUT', 'ASPT', 'ASPT', 'ASPT', 'GLUT', 'GLUT', 'ASPT', 'GLUT', 'ASPT', 'GLUT', 'HSPT', 'GLUT', 'GLUT', 'HSPT', 'GLUT']
    acidics.append('BUF')
    print(acidics)  # debug

    lidx = 0
    for atomset in range(0, len(acidics)):

        if acidics[atomset] == 'HSPT':
            newCoord = str(values[lidx]) + ' ' + str(values[lidx + 1]) + ' ' + str(values[lidx + 2])
            lidx += 3

        if acidics[atomset] in ['ASPT', 'GLUT', 'BUF']:
            newCoord = str(values[lidx])
            lidx += 1

        print(acidics[atomset], newCoord)  # debug

        match   = 'lambda-dynamics-atom-set{}-initial-lambda'.format(atomset + 1)
        replace = match + '               = {}'.format(newCoord)

        # print(match)    # debug
        # print(replace)  # debug

        # Replace the line that contains A with B:
        os.system(f'sed -i \'s/.*{match}.*/{replace}/\' {new}.mdp')


# EM (0.01 * 10000 = 100ps) (calibration = yes)
gromacs('grompp -f EM.mdp -c phneutral.pdb -p topol.top -n index.ndx -o EM.tpr -maxwarn 1')
gromacs(f'mdrun -v -deffnm EM -nt {cores} -npme 0')

# NVT (0.0002 * 50000 = 10ps) (calibration = no)
gromacs('grompp -f NVT.mdp -c EM.gro -p topol.top -n index.ndx -o NVT.tpr -r EM.gro -maxwarn 1')
gromacs(f'mdrun -v -deffnm NVT -nt {cores} -npme 0 -bonded gpu')
extractAndUpdate('NVT', 'HEAVY')

# HEAVY (0.002 * 500000 = 1ns) (calibration = no)
gromacs('grompp -f HEAVY.mdp -c NVT.gro -p topol.top -n index.ndx -o HEAVY.tpr -r NVT.gro -maxwarn 1')
gromacs(f'mdrun -v -deffnm HEAVY -nt {cores} -npme 0 -bonded gpu')
extractAndUpdate('HEAVY', 'BACK')

# BACKBONE (0.002 * 5000000 = 10ns) (calibration = no)
gromacs('grompp -f BACK.mdp -c HEAVY.gro -p topol.top -n index.ndx -o BACK.tpr -r HEAVY.gro -maxwarn 1')
gromacs(f'mdrun -v -deffnm BACK -nt {cores} -npme 0 -bonded gpu')
extractAndUpdate('BACK', 'CA')

# CA (0.002 * 5000000 = 10ns) (calibration = no)
gromacs('grompp -f CA.mdp -c BACK.gro -p topol.top -n index.ndx -o CA.tpr -r BACK.gro -maxwarn 1')
gromacs(f'mdrun -v -deffnm CA -nt {cores} -npme 0 -bonded gpu')
extractAndUpdate('CA', 'MD')
