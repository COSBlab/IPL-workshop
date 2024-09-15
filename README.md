# IPL-workshop
Material for the Structural Bioinformatics workshop at Institut Pasteur Laos (Sept. 30 - Oct. 4 2024).

## Software required
Here is a list of software to be installed on your Linux/Mac laptop.
Some software require a (free) license:

* [Anaconda](https://www.anaconda.com/): management python modules
* [Python Notebook](https://jupyter.org): interactive environment to run python code
* [ChimeraX](https://www.cgl.ucsf.edu/chimerax/): structure visualization
* [FoldX](https://foldxsuite.crg.eu/): scoring function for protein-protein interaction and more
* [DeepRank](https://pypi.org/project/deeprank/): scoring function for protein-protein interaction and more 
* [MODELLER](https://salilab.org/modeller/): homology modelling

### Install python libraries with Conda
Conda can be used to install the python libraries that will be used during the workshop.
First, we need to create a conda environment:
```
conda create --name IPL-workshop
```
and activate it:
```
conda activate IPL-workshop
```
Now we can install all the modules needed during our workshop
```
conda install matplotlib jupyterlab numpy biopython mdanalysis --channel conda-forge
```

Many of the python scripts that we will use will be executed using Jupyter Python Notebooks.
To open the Python Notebooks environment, just type:
```
jupyter lab
```

## Online resources
These are the webservers and online resources that we will use during the workshop:
* [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi): find similarities between sequences 
* [MAFFT](https://mafft.cbrc.jp/alignment/server/index.html): sequence alignment tool
* [AlphaFold2](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb): structure prediction from sequence
* [AlphaFold3](https://alphafoldserver.com/about): structure prediction from sequence
* [Foldseek](https://search.foldseek.com/search): find similarities from structure
* [Pesto](https://pesto.epfl.ch/): parameter-free geometric deep learning method to predict protein interaction interfaces
* [ScanNet](http://bioinfo3d.cs.tau.ac.il/ScanNet/): geometric deep learning model for predicting protein binding sites from structure 
* [HADDOCK](https://www.bonvinlab.org/education/HADDOCK24/HADDOCK24-antibody-antigen-basic/): High Ambiguity Driven protein-protein DOCKing
* [SAbDab](https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab): database containing all the antibody structures available in the PDB

## Program of the workshop
**L**: lecture; **P**: practical

### [Day 1](DAY-1/README.md)
9.00-10.30
*   **L1**:  Introduction to the course: goals and overview of the program
*   **L2**:  Identify homolog structures in complex with receptors and antibodies
*   **L3**:  Revise sequence alignment

10.30-12.30 
*   **P1**: [Identify homolog structures wth BLAST](DAY-1/README.md#P1)
*   **P2**: [Sequence alignment with MAFFT](DAY-1/README.md#P2)

14.00-15.30
*   **L4**: Introduction to MODELLER
*   **L5**: Structural analysis with MDAnalysis

15.30-18.00
*   **P3**: [Build homology models of viral protein in complex with receptor/antibody](DAY-1/README.md#P3)
*   **P4**: [Structural analysis with MDanalysis](DAY-1/README.md#P4)

### [Day 2](DAY-2/README.md) 
9.00-10.00 
*   **L6**: Structure prediction with AI: AlphaFold2 and AlphaFold3
*   **L7**: Foldseek: how to look for similar structures in the PDB
 
10.00-12.30
*   **P5**: [Building models of viral proteins with AF2/AF3](DAY-2/README.md#P5)
*   **P6**: [Look for similar structures in complex with antibodies/receptors](DAY-2/README.md#P6)

14.00-15.00
*   **L8**: Building models of protein complexes with AF3
*   **L9**: Scoring functions to evaluate quality of protein complexes and binding affinity prediction

15.00-18.00
*   **P7**: [Building models of protein complexes with AF3](DAY-2/README.md#P7) 
*   **P8**: [Quality assessment with different scoring functions](DAY-2/README.md#P8)
*   **P9**: [Structural analysis with MDAnalysis](DAY-2/README.md#P9)

### [Day 3](DAY-3/README.md)
9.00-10.00 
*    **L10**: AI prediction of protein binding sites: Pesto and ScanNet
*    **L11**: SAbDab: the Oxford database of antibodies structures

10.00-12.30
*    **P10**: [Practice on the use of Pesto and ScanNet](DAY-3/README.md#P10)

14.00-15.00
*    **L12**: Introduction to physics-based docking: HADDOCK
*    **L13**: Description of the tutorial

15.00-18.00
*   **P11**: [HADDOCK Tutorial on antibody/antigen](DAY-3/README.md#P11)

### [Day 4](DAY-4/README.md)
All day
   * [Group activity: explain goal, provide sequence HEV capsid, go!](DAY-4/README.md#group)
 
### [Day 5](DAY-5/README.md) 
9.00-12.00
* Prepare slides and present results (short group presentation)

