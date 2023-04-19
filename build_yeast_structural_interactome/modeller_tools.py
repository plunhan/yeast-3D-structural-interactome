import pandas as pd
from subprocess import call
from multiprocessing import Process
from modeller import *
from modeller.automodel import *

def produce_protein_models (inPath,
                            alignmentDir,
                            templateDir,
                            modelDir,
                            numModels = 1,
                            verbosity = 'minimal',
                            modellerTimeout = None):
    
    templateMap = pd.read_table (inPath, sep='\t')
    n = len(templateMap)
    
    for i, row in templateMap.iterrows():
        print('\n********************************************************************')
        print('Protein complex: %s, %d out of %d (%.2f%%)' % (row.Complex_ID, i+1, n, 100.*(i+1)/n))
        print('Template file: %s' % row.Template_file_ID)
        print('Alignment file %s' % row.Alignment_file_ID)
        print('********************************************************************\n')
        modelFile = modelDir / (row.Complex_ID + '.B99990001.pdb')
        if not modelFile.is_file():
            delete_model_files (modelDir, row.Complex_ID)
            p = Process (target = create_protein_model,
                         name = "create_protein_model",
                         args = (row.Complex_ID,
                                 row.Template_file_ID,
                                 str(alignmentDir / row.Alignment_file_ID),
                                 str(templateDir),
                                 str(modelDir)),
                         kwargs = {"starting_model":1,
                                   "ending_model":numModels,
                                   "verbosity":verbosity})
            p.daemon = True
            p.start()
            p.join(modellerTimeout)
            if p.is_alive():
                print('Modelling timed out. Killing process...')
                p.terminate()
                p.join()
                delete_model_files (modelDir, row.Complex_ID)

def create_protein_model (protein,
                          templateIDs,
                          alignmentFile,
                          templateDir,
                          modelDir,
                          starting_model = 1,
                          ending_model = 1,
                          verbosity = 'minimal'):
    
    if verbosity is 'verbose':
        log.verbose()
    elif verbosity is 'minimal':
        log.minimal()
    else:
        log.none()
    
    env = environ()
    
    # read heteroatoms
    env.io.hetatm = True
    
    # directories for input atom files
    env.io.atom_files_directory = [modelDir, templateDir]

    # Be sure to use 'MyModel' rather than 'automodel' here!
    a = MyModel(env,
                alnfile  = alignmentFile,
                knowns   = templateIDs,
                sequence = protein)

    a.starting_model = starting_model
    a.ending_model = ending_model
    a.make()

# Override MyModel methods
class MyModel (automodel):
    
    def special_patches (self, aln):
        # number of chains in model
        numChains = len(self.chains)
        
        # Rename both chains and renumber the residues in each
        self.rename_segments (segment_ids = [c.name for c in self.chains],
                              renumber_residues = [1] * numChains)
        
        # Rename chain for single chain model:
        if numChains == 1:
            self.chains[0].name = 'A'

def delete_model_files (modelDir, modelID):
    
    filenames = [modelID + '.B99990001.pdb',
                 modelID + '.D00000001',
                 modelID + '.V99990001',
                 modelID + '.ini',
                 modelID + '.rsr',
                 modelID + '.sch']
                 
    for name in filenames:
        filepath = modelDir / name
        if filepath.is_file():
            call(['rm' , str(filepath)])
