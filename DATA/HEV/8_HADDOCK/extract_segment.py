#!/usr/bin/env python
import sys
import numpy as np
import MDAnalysis as mda

# Usage:
# ./extract_segment.py filename chain res0 resN
# Extract a fragment of a given chain between residues res0 and resN.

filename = sys.argv[1]
chain = sys.argv[2]
res0, resN = sys.argv[3], sys.argv[4]

u = mda.Universe(filename)
atoms = u.select_atoms("chainID {} resid {}:{}".format(chain, res0, resN))
atoms.write(filename.rstrip(".pdb")+"_{}-{}.pdb".format(res0, resN), resindex=False)