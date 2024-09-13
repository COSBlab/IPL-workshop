# IPL-workshop
Material for the Structural Bioinformatics workshop at Institut Pasteur Laos (Sept. 30 - Oct. 4 2024).

## Software required
Here is a list of software to be installed on your Linux/Mac laptop.
Some software require a (free) license:

* [Anaconda](https://www.anaconda.com/): management python modules
* [Python Notebook](https://jupyter.org): interactive environment to run python code
* [MDAnalysis](https://www.mdanalysis.org/): python library for structural analysis
* [ChimeraX](https://www.cgl.ucsf.edu/chimerax/): structure visualization
* [FoldX](https://foldxsuite.crg.eu/): scoring function for protein-protein interaction and more
* [DeepRank](https://pypi.org/project/deeprank/): scoring function for protein-protein interaction and more 
* [MODELLER](https://salilab.org/modeller/): homology modelling

## Online resources
These are the webservers and online resources we will use during the workshop:
* [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi): find similarities between sequences 
* [MAFFT](https://mafft.cbrc.jp/alignment/server/index.html): sequence alignment tool
* [AlphaFold2](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb): structure prediction from sequence
* [AlphaFold3](https://alphafoldserver.com/about): structure prediction from sequence
* [FoldSeek](https://search.foldseek.com/search): find similarities from structure
* [Pesto](https://pesto.epfl.ch/): parameter-free geometric deep learning method to predict protein interaction interfaces
* [ScanNet](http://bioinfo3d.cs.tau.ac.il/ScanNet/): geometric deep learning model for predicting protein binding sites from structure 
* [HADDOCK](https://www.bonvinlab.org/education/HADDOCK24/HADDOCK24-antibody-antigen-basic/): High Ambiguity Driven protein-protein DOCKing
* [SAbDab](https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab): database containing all the antibody structures available in the PDB

## Preliminary program
**T**: theoretical lecture; **P**: practical

### Day 1
9.00-10.30
*   **T**:  Introduction to the course: goals and overview of the program
*   **T**:  How to BLAST the PDB for structures with partners
*   **T**:  Revise sequence alignment (with ClustalW or MAFFT)

10.30-12.30 
*   **P**: Find homolog structures with antibodies/receptors + sequence alignment: Exercise. Test system: XXX

14.00-15.30
*   **T**: Introduction to MODELLER
*   **T**: Structural analysis with MDAnalysis

15.30-18.00
*   **P**: Build homology models viral protein + receptor/antibody. Test system: XXX
*   **P**: Structural analysis (map interface) with MDanalysis (We can provide a jupyter notebook template for interface mapping to use each time they need to do this if they have no prior notion of python)

### Day 2
9.00-10.00 
*   **T**:  Introduction to AF2 and AF3 for monomer modeling
*   **T**:  FoldSeek: how to look for similar structures in the PDB (with receptors/antibodies) and retrieve data
 
10.00-12.30
*   **P**: Building models of viral proteins with AF2/AF3. Test system: XXX
*   **P**: Look for similar structures in complex with antibodies/receptors
*   **P**: Structural analysis (map interface) with MDAnalysis

14.00-15.00
*   **T**: Building models of protein complexes with AF2/AF3-multimer
*   **T**: Scoring functions to evaluate quality of protein complexes and binding affinity prediction

15.00-18.00
*   **P**: AF2/AF3-multimer (viral protein + partners identified before from homology or FoldSeek). Test system: XXX
*   **P**: Structural analysis (map interface) with MDAnalysis
*   **P**: Quality assessment with different scoring functions

### Day 3
9.00-10.00 
*    **T**: AI prediction of protein binding sites: Pesto and ScanNet
*    **T**: The Oxford database of antibodies structures

10.00-12.30
*    **P**: Practice on the use of Pesto and ScanNet. Test Systems: XXX
*    **P**: Benchmark on bound structures of antigen/antibodies from Oxford database

14.00-15.00
*    **T**:  Introduction to physics-based docking: HADDOCK
*    **T**:  Description of the tutorials

15.00-18.00
*   **P**: HADDOCK Tutorial antibody/antigen

### Day 4: group activity
All day
   * Explain goal, provide data
   * Start exercise
 
### Day 5
9.00-12.00
* prepare slides and present results (short group presentation)

