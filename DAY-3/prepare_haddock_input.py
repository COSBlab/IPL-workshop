#!/usr/bin/env python
import sys
import numpy as np
import MDAnalysis as mda
# example: AF model with chain A antigen and chain B and C entibody
# antibody extraction: ./prepare_haddock_input.py filename B C
# antigen extraction: ./prepare_haddock_input.py filename A


filename = sys.argv[1]
chains = sys.argv[2:]


u = mda.Universe(filename)


#output_atoms
for i, chain in enumerate(chains):
    print(chain)
    if i == 0:
        output_atoms = u.select_atoms("chainID {}".format(chain))
    else:
        max_res = np.max(output_atoms.resids)
        atoms = u.select_atoms("chainID {}".format(chain))
        for j, el in enumerate(atoms):
            el.chainID = chains[0]
            if el.name == "CA":
                el.residue.resid += max_res
        output_atoms += atoms


# save
chains_str = "_"
for i, el in enumerate(chains):
    if i != 0:
        chains_str += "_"
    chains_str += el
outfile = filename.rstrip(".pdb")+chains_str+".pdb"
output_atoms.write(outfile, reindex=False)