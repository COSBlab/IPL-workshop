# IPL-workshop: Day 3

## <a name="P10"></a>P10: AI prediction of protein-protein interactions: PeSTo and ScanNet

**Aim**

In this exercise, we will use [PeSTo](https://pesto.epfl.ch/) and [ScanNet](http://bioinfo3d.cs.tau.ac.il/ScanNet/) to predict protein-protein and protein-antibodies interfaces starting from a structural model of the viral protein.  

**Tasks**

To complete this exercise,

* we will first start with a quick benchmark of the two approaches. To do so, we will use known experimental structures of antibody/antigen complexes retrieved from [SAbDab](https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab), a database containing all the antibody structures available in the PDB. Structures satisfying certain criteria have already beem downloaded locally: XXX

* the students need to extract from a selection of experimental structures the structure of the antigen. To do so, the student can use the ```extract_chains.py``` script utilized in one of our previous exercises. The table ```XXX``` contains a list of PDBs along with chains assignement for antibody and antigen.

* The structures of the antigens will be uploaded to PeSTo and ScanNet servers:
  1. for PeSTo: choose the model ```PeSTo```, select ```Detect chains```, and then ```Submit```
  2. for ScanNet: choose binding site type ```Protein-antibody```,  provide an email address, and then ```Submit```

* In both cases, the students will retrieve a PDB file in which the Bfactor column contains a scoring function quantifying the probability of each residue to participate in protein-protein (or protein-antibody) interactions. To analyze the PDBs, the student will:
  1. visual the PDBs with ChimeraX and color each residue by Bfactor 
  2. adapt our previous Python Notebooks to create a plot showing, for each residue, the value of the PeSTo and ScanNet scoring function and comparing these predictions with the interacting residues extracted from the PDB of the antibody/antigen complex deposited in [SAbDab](https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab)

* After benchmarking the two methods on a few complexes extracted from SAbDab, the student will make predictions for our viral protein starting from their best AF2/AF3 models

* As for the SAbDab benchmark, the student will adapt our previous Python Notebooks to create a plot showing, for each residue, the value of the PeSTo and ScanNet scoring function and comparing these predictions with the interacting residues in models of viral protein/antibodies complexes obtained with MODELLER and AF3 in previous exercises.


## <a name="P11"></a>P11: HADDOCK Tutorial on antibody/antigen

**Aim**

In this exercise, we will learn how to use [HADDOCK](https://rascar.science.uu.nl/haddock2.4/) to model antibody/antigen complexes. 

**Tasks**

To complete this exercise, the student will:

* complete the tutorial on antibody/antigen available [here](https://www.bonvinlab.org/education/HADDOCK24/HADDOCK24-antibody-antigen-basic/). Credentials to access the server will be distributed during the workshop
