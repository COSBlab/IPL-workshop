# IPL-workshop: Day 2

## P5: Building models of viral proteins with AF2/AF3

**Aim**

In this exercise, we will use both [AF2](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb) and [AF3](https://alphafoldserver.com/about) to generate structural models of the viral protein directly from sequence and analyze the quality of the models.

**Tasks**

To complete this exercisei:
* we will start using the ColabFold v1.5.5 implementation of AF2 to generate structural models.
  The students is expected to play with the following parameters:
  1. ```template_mode```: to avoid using templates (```none```) or to use detected templates in pdb100 (```pdb100```)
  2. ```num_recycles```: to change the number of recycles (try varying from 1 to 3)
  3. playing with sampling parameters: ```max_msa```, ```num_seeds```, and ```use_dropout``` 

* The quality of the models (avg.pLDDT, MSA coverage) can be inspected directly online by inspecting the plots reported. 

  We will download the PDB of the models for visualization with ChimeraX and analysis with Python notebook (see below)

* we will then use AF3, which provides a limited number of options. With AF3, we will inspect the model quality online and download a zip file with all the results for further analysis.

  ![title](AF3.png)

* now that we have generated a number of models with AF2 and AF3, we will perform additional analysis on our local hardware.
  1. we start with visual inspection. You can open all the models with ChimeraX, align then with "Tools/Structure Analysis/Matchmaker" and color them by Bfactors. We can use the AlphaFold classic palette, by typing
     ```color bfactor palette alphafold```` in ChimeraX command line. 

## P6: Look for similar structures in complex with antibodies/receptors

**Aim**

**Tasks**

