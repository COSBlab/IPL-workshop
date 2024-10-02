# IPL-workshop

Material for the Structural Bioinformatics workshop at Institut Pasteur Laos (Sept. 30 - Oct. 4 2024).

Authors: Vincent Schnapka and Max Bonomi (Institut Pasteur - CNRS).

## Software required

Here is a list of software to be installed on your Linux/Mac laptop:

* [Anaconda](https://www.anaconda.com/): management python modules
* [FoldX 5.1](https://foldxsuite.crg.eu/): scoring function for protein-protein interaction and more
* [ChimeraX](https://www.cgl.ucsf.edu/chimerax/): structure visualization
* [PyMOL](https://www.pymol.org/): structure visualization (optional)

**Note** Some software require a (free) academic license. Please make sure you obtain one before joining the workshop.

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
And add MODELLER:
```
conda install modeller --channel salilab
```

Many of the python scripts that we will use will be executed using Jupyter Python Notebooks.
To open the Python Notebooks environment, just type:
```
jupyter lab
```

**Note** Please make sure to setup the conda environment and install the libraries prior to the meeting.

## Online resources

These are the webservers and online resources that we will use during the workshop:

* [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi): find similarities between sequences 
* [MAFFT](https://mafft.cbrc.jp/alignment/server/index.html): sequence alignment tool
* [AlphaFold2](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb): structure prediction from sequence
* [AlphaFold3](https://alphafoldserver.com/about): structure prediction from sequence
* [Foldseek](https://search.foldseek.com/search): find similarities from structure
* [PeSTo](https://pesto.epfl.ch/): parameter-free geometric deep learning method to predict protein interaction interfaces
* [ScanNet](http://bioinfo3d.cs.tau.ac.il/ScanNet/): geometric deep learning model for predicting protein binding sites from structure 
* [HADDOCK](https://www.bonvinlab.org/education/HADDOCK24/HADDOCK24-antibody-antigen-basic/): High Ambiguity Driven protein-protein DOCKing
* [SAbDab](https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab): database containing all the antibody structures available in the PDB

**Note** The AlphaFold3 server requires opening an account. Please make sure you set this up before joining the workshop.
For HADDOCK, we will provide credentials to use the server. However, please consider opening your own account prior to the meeting.

## Program of the workshop

**L**: lecture; **P**: practical

**Coffee breaks: 10.00-10.30 and 15.00-15.30**

### [Day 1](DAY-1/README.md)

8.30-12.00

*   **L1**:  Introduction to the course: goals and overview of the program
*   **L2**:  Identify homolog structures in complex with receptors and antibodies
*   **P1**: [Identify homolog structures wth BLAST](DAY-1/README.md#P1)
*   **L3**:  Revise sequence alignment
*   **P2**: [Sequence alignment with MAFFT](DAY-1/README.md#P2)

12.00-13.00 Lunch

13.00-16.30

*   **L4**: Introduction to MODELLER
*   **P3**: [Build homology models of viral protein in complex with receptor/antibody](DAY-1/README.md#P3)
*   **L5**: Structural analysis with MDAnalysis
*   **P4**: [Structural analysis with MDAnalysis](DAY-1/README.md#P4)

16.30-17.00 Debriefing and discussions

### [Day 2](DAY-2/README.md) 

8.30-12.00

*   **L6**: Structure prediction with AI: AlphaFold2 and AlphaFold3
*   **P5**: [Building models of viral proteins with AF2/AF3](DAY-2/README.md#P5)
*   **L7**: Foldseek: how to look for similar structures in the PDB
*   **P6**: [Look for similar structures in complex with antibodies/receptors](DAY-2/README.md#P6)

12.00-13.00 Lunch

13.00-16.30

*   **L8**: Structure prediction of protein complexes with AlphaFold
*   **P7**: [Building models of protein complexes with AF3](DAY-2/README.md#P7)
*   **L9**: Binding free-energy estimation for protein complexes
*   **P8**: [Quality assessment and structural analysis with MDAnalysis](DAY-2/README.md#P8)
*   **P9**: [Binding free-energy estimation](DAY-2/README.md#P9)

16.30-17.00 Debriefing and discussions

### [Day 3](DAY-3/README.md)

8.30-12.00

*   **L10**: Predicting protein-protein interactions with AI
*   **L11**: SAbDab: the Oxford database of antibodies structures
*   **P10**: [Predicting protein-protein interactions with AI](DAY-3/README.md#P10)
*   **L12**: Introduction to physics-based docking: HADDOCK

12.00-13.00 Lunch

13.00-15.30

*   **P11**: [HADDOCK Tutorial on antibody/antigen](DAY-3/README.md#P11)

15.30-15.50 Debriefing and discussions

### [Day 4](DAY-4/README.md)

8.30-16.30

* [Group activity: explain goal, provide sequence HEV capsid, go!](DAY-4/README.md#group)
 
16.30-17.00 Debriefing and discussions

### [Day 5](DAY-5/README.md) 

8.30-11.30

* [Prepare slides and present results](DAY-5/README.md#results)
* Identify constructs for HEV LIPS
* Feedback on the workshop

11.30-12.00 Attendance certificates
