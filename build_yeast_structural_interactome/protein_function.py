#----------------------------------------------------------------------------------------
# Modules for computations on protein function.
#----------------------------------------------------------------------------------------

import pickle
import pandas as pd
import subprocess

def produce_fastsemsim_protein_gosim_dict (inPath,
                                           outPath,
                                           task = 'SS',
                                           ont_type = 'GeneOntology',
                                           sim_measure = 'SimGIC',
                                           mix_method = 'BMA',
                                           cutoff = -1,
                                           remove_nan = True,
                                           query_mode = 'pairs',
                                           query_ss_type = 'obj',
                                           ac_species = 'human',
                                           ont_root = 'biological_process',
                                           ont_ignore = None,
                                           ec_ignore = None,
                                           verbosity = '-vv',
                                           ontologyFile = 'go-basic.obo',
                                           annotationFile = None,
                                           paramOutFile = 'fastsemsim_parameters',
                                           fastsemsimOutFile = 'fastsemsim_output'):
    """Calculate gene ontology (GO) similarity using fastsemsim library.

    Args:
        
    
    Returns:
        
    
    """
    if annotationFile:
        cmd = ['fastsemsim', '--ac_file', annotationFile]
    else:
        cmd = ['fastsemsim', '--ac_species', ac_species]
    cmd += [verbosity, '--task', task, '-o', ont_type, '--ontology_file', ontologyFile, 
            '--query_input', 'file', '--query_mode', query_mode, '--query_ss_type', query_ss_type, 
            '--query_file', inPath, '--tss', sim_measure, '--tmix', mix_method, '--root', 
            ont_root, '--cut', str(cutoff), '--remove_nan', '--save_params', paramOutFile, 
            '--output_file', fastsemsimOutFile]
    if ec_ignore:
        for ec in ec_ignore:
            cmd += ['--ignore_EC', ec]
    if ont_ignore:
        for ont in ont_ignore:
            cmd += ['--ontology_ignore', ont]
    if remove_nan:
        cmd.append('--remove_nan')
    print(' '.join(map(str, cmd)))
    subprocess.run(cmd)
    table = pd.read_table(fastsemsimOutFile, sep='\t')
    gosim = {tuple(sorted((p1, p2))):sim for p1, p2, sim in table.values}
    with open(outPath, 'wb') as fOut:
        pickle.dump(gosim, fOut)
