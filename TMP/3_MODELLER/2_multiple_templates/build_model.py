#!/usr/bin/env python
from modeller import *
from modeller.automodel import *


env4 = Environ()
env4.io.atom_files_directory = ['.', '../../Structures/']

a = AutoModel(env4, 
              alnfile='alignment.pir',
              knowns=[''], 
              sequence='capsid',
              assess_methods=(assess.DOPE,
                              assess.normalized_dope,
                              assess.GA341))
a.starting_model = 1
a.ending_model = 5
a.make()
