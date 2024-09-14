# Multiple Sequence Alignment with MAFFT

provided:
- A jupyter notebook to make the input fasta for MAFFT (make_input.ipynb)
- a pdb2fasta python conversion script (pdb2fasta.py)

provide:
A BLAST hit table containing the list of homologous PDBs (hit_table.txt)

step1:
In the jupyter notebook: 
- Read the hit table and select the best homologous sequences
- Make the fasta sequences and prepare the input for MAFFT (homolog_sequences.fasta)

step2:
Run MAFFT in the webserver to perform MSA of the homologous sequences

step3:
Download the results in CLUSTAL format and fasta format.
