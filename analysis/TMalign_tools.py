
def retrieve_mapping_data(output): 
    # Input (str): TM-align output information
    # Output (str, str, str): Residues and gaps for sequences 1 and 2, and mapping details
    lines = output.decode("utf-8").split('\n')
    mark = '(":" denotes residue pairs of d <  5.0 Angstrom, "." denotes other aligned residues)'
    ind = lines.index(mark)
    seq_1 = lines[ind+1]
    seq_2 = lines[ind+3]
    mapping = lines[ind+2]
    return seq_1, seq_2, mapping

def assign_residue_position(seq): 
    # Input (str): Residues and gaps for a sequence
    # Output (list): A list of positions of residues
    i = 0
    outputList = []
    for res in seq: 
        if res != '-': 
            i += 1
            outputList.append(i)
        else: 
            outputList.append('-')
    return outputList

def check_mapping_status_for_residue(resPos, outputList, mapping): 
    # Input: 
    # resPos (int): The position of residue to be checked
    # outputList (list): The list of positions of residues
    # mapping (list): The list of characters showing mapping details
    # Output (boolean): mapped (True) or not (False)
    ind = outputList.index(resPos)
    if mapping[ind] in [':', '.']: 
        return True
    else: 
        return False

def check_mapping_status_for_interface(output, intResPos_1, intResPos_2, strSim_cutoff): 
    seq_1, seq_2, mapping = retrieve_mapping_data(output)
    resPos_1 = assign_residue_position(seq_1)
    resPos_2 = assign_residue_position(seq_2)
    mapping_1, mapping_2 = 0, 0
    for intResPos in intResPos_1: 
        if check_mapping_status_for_residue(intResPos, resPos_1, mapping): 
            mapping_1 += 1
    for intResPos in intResPos_2: 
        if check_mapping_status_for_residue(intResPos, resPos_2, mapping): 
            mapping_2 += 1
    if mapping_1 / len(intResPos_1) >= strSim_cutoff and mapping_2 / len(intResPos_2) >= strSim_cutoff: 
        return True
    else: 
        return False