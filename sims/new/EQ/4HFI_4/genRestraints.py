from science.utility import genRestraints

# Note: we don't have to generate HEAVY.itp as this is simply the original posre.itp from pdb2gmx.
genRestraints('phneutral.pdb', 'BACK.itp', 'chainID A and not resname POPC and name CA C N O')  # good
genRestraints('phneutral.pdb', 'CA.itp',   'chainID A and not resname POPC and name CA')        # good
