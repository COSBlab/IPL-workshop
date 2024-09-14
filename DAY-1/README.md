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

* retrieve the sequences of the homologs identified in the previous exercise with BLAST
* paste these sequences along with the target sequence in the "Input" box
* run the alignment with default parameters
* inspect the results
* download the alignment file (Fasta format) to be used in the next exercises 
