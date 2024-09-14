import sys
# Homology modeling with multiple templates
from modeller import *              # Load standard Modeller classes
from modeller.automodel import *    # Load the automodel class

# create enviroment
log.verbose()    # request verbose output
env = environ()  # create a new MODELLER environment to build this model in

# directories for input atom files
env.io.atom_files_directory = ['.', './']

# define parameters
a = AutoModel(env,
              alnfile  = 'alignment.pir', # alignment filename
              knowns   = ('template'),  # codes of the templates
              sequence = 'BANAL-20-52')     # code of the target
# number of models
a.starting_model= 1  # index of the first model
a.ending_model  = 5  # index of the last model
# scoring functions
a.assess_methods = (assess.DOPE,
                    assess.normalized_dope)
# do the actual homology modeling
a.make()
