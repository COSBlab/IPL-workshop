import MDAnalysis as mda
import sys

# PDB FILE
FILE_=sys.argv[1]
# chains to extract
ch = sys.argv[2:]

# load universe template PDB
u = mda.Universe(FILE_)

# get atoms of desired chains
at_p = u.select_atoms("segid "+" ".join(ch))

# write a PDB with desired chains 
at_p.write("template.pdb")
