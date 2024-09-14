import MDAnalysis as mda
import sys

# PDB FILE
FILE_=sys.argv[1]

# dictionary 3-to-1 code
al={}
al["ALA"]="A"; al["ARG"]="R"; al["ASN"]="N"; al["ASP"]="D"; al["CYS"]="C"
al["GLN"]="Q"; al["GLU"]="E"; al["GLY"]="G"; al["HSD"]="H"; al["ILE"]="I"
al["LEU"]="L"; al["LYS"]="K"; al["MET"]="M"; al["PHE"]="F"; al["PRO"]="P"
al["SER"]="S"; al["THR"]="T"; al["TRP"]="W"; al["TYR"]="Y"; al["VAL"]="V"
al["HIS"]="H"; al["HIE"]="H"

# load universe reference PDB complex
u = mda.Universe(FILE_)

# get protein atoms
at_p = u.select_atoms("protein")
# get list of chains
ch = at_p.segments

# prepare output file
out = open(FILE_.split(".")[0]+".fasta", "w")

# loop over chains
for c in ch:
    # initialize sequence
    seq=''
    # loop over residues
    for r in at_p.select_atoms('segid '+str(c.segid)).residues:
        seq += al[r.resname]
    out.write(">chain-%s\n" % str(c.segid))
    out.write("%s\n" % seq)

# write a PDB with only protein atoms
u.select_atoms("protein").write(FILE_.split(".")[0]+"-protein.pdb")
