# IPL-workshop: Day 2

## P5: Building models of viral proteins with AF2/AF3

**Aim**

In this exercise, we will use both [AF2](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb) and [AF3](https://alphafoldserver.com/about) to generate structural models of the viral protein directly from sequence and analyze the quality of the models.

**Tasks**

To complete this exercise:
* we will start using the ColabFold v1.5.5 implementation of AF2 to generate structural models.
  The students might want to play with the following parameters:
  1. ```template_mode```: to avoid using templates (```none```) or to use detected templates in pdb100 (```pdb100```)
  2. ```num_recycles```: to change the number of recycles (try varying from 1 to 3)
  3. the sampling parameters: ```max_msa```, ```num_seeds```, and ```use_dropout``` 

* The quality of the models (Predicted aligned error, or PAE, MSA coverage, and pLDDT) can be assessed directly online by inspecting the plots reported. 

  ![title](AF2-plots.png)

  We will download the PDB of the models for visualization with ChimeraX and analysis with Python notebook (see below)

* we will then use AF3, which does not provide the same options as the ColabFold implementation of AF2. With AF3, we will inspect the model quality online and download a zip file with all the results for further analysis.

  ![title](AF3.png)

* now that we have generated a number of models with AF2 and AF3, we will perform additional analysis on our local hardware.
  1. we start with visual inspection. You can open all the models with ChimeraX, align then with "Tools/Structure Analysis/Matchmaker" and color them by Bfactors. We can use the AlphaFold classic palette, by typing ```color bfactor palette alphafold``` in ChimeraX command line. 

  2. now we will analyse the per-residue pLDDT across all models generated with AF2 and AF3. To do so, we have prepared a Python Notebook:

     ```jupyter lab analyze_AF.ipynb```

  3. based on the quality assessment above, we will identify the most reliable model(s) to be used in the next exercise 

## P6: Look for similar structures in complex with antibodies/receptors

**Aim**

In this exercise, we will use [FoldSeek](https://search.foldseek.com/search) to identify proteins structurally similar to the AF2/AF3 of the viral spike and possibly in complex with antibodies/receptors. Sometimes such proteins have low sequence identity with the target protein and therefore cannot be identified with a sequence similarity search with BLAST.

**Tasks**

To complete this exercise, the student will:
* select one or more AF2/AF3 models obtained in the previous exercise
* upload the PDB of the model on FoldSeek. For AF3 models, please make sure to use the PDB file and not the CIF.
* select the ```PDB100`` database
* try with two different similarity metrics: ```3Di/AA``` and ```TM-align```
* inspect the table of results, in particular: 

