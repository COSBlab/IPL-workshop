# IPL-workshop: Day 3

## <a name="P10"></a>P10: Predicting protein-protein interactions with AI: PeSTo and ScanNet

**Aim**

In this exercise, we will use [PeSTo](https://pesto.epfl.ch/) and [ScanNet](http://bioinfo3d.cs.tau.ac.il/ScanNet/) to predict protein-protein and protein-antibodies interfaces starting from a structural model of the viral protein.  

**Tasks**

To complete this exercise,

* we will first start with a quick benchmark of the two approaches. To do so, we will use known experimental structures of antibody/antigen complexes retrieved from [SAbDab](https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab), a database containing all the antibody structures available in the PDB. 52 structures satisfying certain criteria have already beem downloaded locally: protein antigen, max sequence identity 65%, maximum resolution 3 Ang.

* the students need to extract from a selection of experimental structures (around 5-10) the structure of the antigen. To do so, the student can use the ```extract_chains.py``` script utilized in one of our previous exercises. The table ```summary.tsv``` contains a list of PDBs along with chains assignement for antibody and antigen.

* The structures of the antigens will be uploaded to PeSTo and ScanNet servers:
  1. for PeSTo: choose the model ```PeSTo```, select ```Detect chains```, and then ```Submit```
  2. for ScanNet: choose binding site type ```Protein-antibody```,  provide an email address, and then ```Submit```

* In both cases, the students will retrieve a PDB file in which the Bfactor column contains a scoring function quantifying the probability of each residue to participate in protein-protein (or protein-antibody) interactions. To analyze the PDBs, the student will:
  1. visualize the PDBs with ChimeraX and color each residue by Bfactor 
  2. use the Notebook ```analyze_Pesto_Scannet.ipynb``` to show, for each residue, the value of the PeSTo and ScanNet scoring function and to compare these predictions with the interacting residues extracted from the PDB of the antibody/antigen complex deposited in [SAbDab](https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab)

* After benchmarking the two methods on a few complexes extracted from SAbDab, the student will make predictions for our viral protein starting from their best AF2/AF3 models.
  If needed, a model in ```cif``` format called ```model.cif``` can be converted into a ```pdb``` file on the command line:
  ```
  python convert_cif_pdb.py model.cif
  ```

* The students will compare these predictions with the results obtained in the previous exercises (visually or using the Notebook ```analyze_Pesto_Scannet.ipynb```)


## <a name="P11"></a>P11: HADDOCK Tutorial on antibody/antigen

**Aim**

In this exercise, we will learn how to use [HADDOCK](https://rascar.science.uu.nl/haddock2.4/) to model antibody/antigen complexes. 

**Tasks**

To complete this exercise, the student will:

* complete the tutorial on antibody/antigen available [here](https://www.bonvinlab.org/education/HADDOCK24/HADDOCK24-antibody-antigen-basic/). Credentials to access the server will be distributed during the workshop.

* In the visualisation part of the tutorial, if PyMOL crashes upon running the rms_cur command, the student can use the provided "calc_lrmsd" notebook to calculate the ligand-rmsd of the different models.
