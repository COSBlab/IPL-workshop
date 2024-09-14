#!/usr/bin/env python
import numpy as np
from argparse import ArgumentParser


__version__ = "0.0.1"


def build_parser():
    parser = ArgumentParser(description="Convert protein pdb file into the corresponding fasta sequence")
    parser.add_argument('-pdb', type=str, help="Input alphafold PDB or CIF file")
    parser.add_argument('-out', type=str, help="Output fasta file (default: protein.fasta)", default="protein.fasta")
    parser.add_argument('-chain', type=str, help="Chain ID (default=A)", default="A")
    return parser


seq = {"G": "GLY", "R": "ARG", "K": "LYS", "D": "ASP", "E": "GLU", "S": "SER", "T": "THR", "N": "ASN", "Q": "GLN", "C": "CYS", "U": "SEC",
       "P": "PRO", "A": "ALA", "V": "VAL", "I": "ILE", "L": "LEU", "M": "MET", "F": "PHE", "Y": "TYR", "W": "TRP", "H": "HIS"}


inv_seq = {'GLY': 'G',
 'ARG': 'R',
 'LYS': 'K',
 'ASP': 'D',
 'GLU': 'E',
 'SER': 'S',
 'THR': 'T',
 'ASN': 'N',
 'GLN': 'Q',
 'CYS': 'C',
 'SEC': 'U',
 'PRO': 'P',
 'ALA': 'A',
 'VAL': 'V',
 'ILE': 'I',
 'LEU': 'L',
 'MET': 'M',
 'PHE': 'F',
 'TYR': 'Y',
 'TRP': 'W',
 'HIS': 'H',
 'HIE': 'H'}


def read_pdb(input="protein.pdb", chain="A"):
    # Extract the data:
    with open(input, "r") as f:
        data = [el.rstrip("\n") for el in f.readlines() if el.split()[0] == 'ATOM']
    # Atom numbers in the PDB file:
    atom_idx = [(int(el[7:11])) for el in data]
    # Atom types in the PDB file:
    atom_type = [el[12:17].rstrip(" ").lstrip(" ") for el in data]
    # Residue type in the PBD file:
    res_type = [el[17:20] for el in data]
    # Residue indexes in the PDB file:
    res_idx = [int(el[23:26]) for el in data]
    # chain id:
    chain_idx = [el[21] for el in data]
    # Storing the atom number of CB (or CA if GLY) and the residue type for each residue:
    res_data = {k: (atom_idx[i], res_type[i]) for i, k in enumerate(res_idx) if chain_idx[i] == chain and (atom_type[i] == "CB" or (atom_type[i] == "CA" and res_type[i] == "GLY"))}
    return res_data


def pdb2sequence(input="protein.pdb", chain="A"):
    # Extract the data:
    with open(input, "r") as f:
        data = [el.rstrip("\n") for el in f.readlines() if el.split()[0] == 'ATOM']
    # Residue type in the PBD file:
    res_idx = [int(el[23:26]) for el in data]
    res_type = [el[17:20] for el in data]
    chain_idx = [el[21] for el in data]
    res_dic = {idx: res_type[i] for i, idx in enumerate(res_idx) if chain_idx[i] == chain}
    sequence = [inv_seq[res_dic[i]] for i in res_dic.keys()]
    return sequence


def pdb2fasta(input_pdb_file="protein.pdb", output_fasta_file="fasta.fasta", chain_ID="A"):
    sequence = pdb2sequence(input=input_pdb_file, chain=chain_ID)
    with open(output_fasta_file, "w") as f:
        pdbname = input_pdb_file.rstrip(".pdb").split("/")[-1]
        f.write(">{}\n".format(pdbname))
        for el in sequence:
            f.write(el)
        #f.write("\n")


def main(pdb_file, output_file, chain_ID):
    pdb2fasta(pdb_file, output_file, chain_ID)


if __name__ == "__main__":
    args = build_parser().parse_args()
    main(args.pdb, args.out, args.chain)
