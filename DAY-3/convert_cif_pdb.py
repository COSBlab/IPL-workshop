from Bio import PDB
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.StructureBuilder import StructureBuilder
import re
import sys

# cif filename
CIF_=sys.argv[1]

# get PDB file name
fpdb=re.sub(".cif",".pdb",CIF_)
# convert to PDB
parser = PDB.MMCIFParser()
structure = parser.get_structure(structure_id="PDB", filename=CIF_)
# save PDB
io=PDBIO()
io.set_structure(structure)
io.save(fpdb)
