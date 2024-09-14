# IPL-workshop: Day 1

## P1: Identify homolog structures wth BLAST

**Aim**

In this exercise, we will use [Protein BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) to look in the PDB database for systems with high sequence identity to the target and possibly in complex with antibodies and receptors.

**Tasks**

To complete this exercise, the student needs to:
* paste the target sequence into the box "Enter Query Sequence"
* select the "Protein Data Bank" as database
* run the search
* analyse the results table: system type, sequence identity, coverage, structure resolution
* identify homologs and inspect them visually with [ChimeraX](https://www.cgl.ucsf.edu/chimerax/).

The results table should look like this:

![title](blast.png)

## P2: Sequence alignment with MAFFT

**Aim**

In this exercise, we will use [MAFFT](https://mafft.cbrc.jp/alignment/server/index.html) to do a sequence alignment between the homolog(s) identified in the previous exercise and the target sequence.

**Tasks**

To complete this exercise, the student needs to:

* retrieve the sequences of the homolog identified in the previous exercise with BLAST. 
  The sequence should be extracted from each PDB file using the ```get_fasta_PDB.py``` script as follows:
```
python get_fasta_PDB.py 7cwl.pdb 
```
This script will create two files:
1. ```7cwl.fasta```: sequence extracted from the PDB, one per chain. Non-protein atoms have been removed
2. ```7cwl-protein.pdb```: PDB file containing only protein residues matching the sequence file ```7cwl.fasta``` 

* paste the correct sequence (ignore sequences of antibodies/receptors present in the homolog PDB) along with the target sequence in the "Input" box
* run the alignment with default parameters
* inspect the results
* download the alignment file (Fasta format) to be used in the next exercises 
