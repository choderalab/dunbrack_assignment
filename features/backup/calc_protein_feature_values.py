"""
protein_features.py
This is a tool to featurize kinase conformational changes through the entire Kinome.

"""

def key_klifs_residues(numbering):
    """
    Retrieve a list of PDB residue indices relevant to key kinase conformations mapped via KLIFS.

    Define indices of the residues relevant to a list of 12 collective variables relevant to
    kinase conformational changes. These variables include: angle between aC and aE helices,
    the key K-E salt bridge, DFG-Phe conformation (two distances), X-DFG-Phi, X-DFG-Psi,
    DFG-Asp-Phi, DFG-Asp-Psi, DFG-Phe-Phi, DFG-Phe-Psi, DFG-Phe-Chi1, and the FRET L-S distance.
    All features are under the current numbering of the structure provided.

    Parameters
    ----------
    numbering : list of int
        numbering[klifs_index] is the residue number for the given PDB file corresponding to KLIFS residue index 'klifs_index'

    Returns
    -------
    key_res : list of int
        Key residue indices

    """
    if numbering == None:
        print("The structure was not found in the klifs database.")
        key_res = None
        return key_res

    key_res = dict() #initialize key_res (which read from the 0-based numbering list)
    for i in range(5):
        key_res[f'group{i}'] = list()
    ## feature group 0: A-loop backbone dihedrals
    key_res['group0'].append(numbering[83]) # start of A-loop

    ## feature group 1: P-loop backbone dihedrals
    key_res['group1'].append(numbering[3]) # res0 in P-loop
    key_res['group1'].append(numbering[4]) # res1 in P-loop
    key_res['group1'].append(numbering[5]) # res2 in P-loop
    key_res['group1'].append(numbering[6]) # res3 in P-loop
    key_res['group1'].append(numbering[7]) # res4 in P-loop
    key_res['group1'].append(numbering[8]) # res5 in P-loop

    ## feature group 2: aC-related features
    #angle between aC and aE helices
    key_res['group2'].append(numbering[19])  # res0 in aC
    key_res['group2'].append(numbering[29])  # res10 in aC
    key_res['group2'].append(numbering[62])  # end of aE

    # key salt bridge
    key_res['group2'].append(numbering[16])  # K in beta III
    key_res['group2'].append(numbering[23])  # E in aC

    ## feature group 3: DFG-related features
    key_res['group3'].append(numbering[79])  # X-DFG
    key_res['group3'].append(numbering[80])  # DFG-Asp
    key_res['group3'].append(numbering[81])  # DFG-Phe
    key_res['group3'].append(numbering[27])  # ExxxX

    ## feature group 4: the FRET distance
    # not in the list of 85 (equivalent to Aura"S284"), use the 100% conserved beta III K as a reference
    key_res['group4'].append(numbering[16] + 120)

    # not in the list of 85 (equivalent to Aura"L225"), use the 100% conserved beta III K as a reference
    key_res['group4'].append(numbering[16] + 61)

    return key_res

def compute_simple_protein_features(traj, key_res):
    """
    This function takes the PDB code, chain id and certain coordinates of a kinase from
    a command line and returns its structural features.

    Parameters
    ----------
    traj : str
	A MDTraj.Trajectory object of the input structure (a pdb file or a simulation trajectory).
    key_res : dict of int
        A dictionary (with keys 'group0' ... 'group4') of feature-related residue indices in five feature groups.
    Returns
    -------
    features: list of floats
    	A list (single structure) or lists (multiple frames in a trajectory) of 72 features in 5 groups (A-loop, P-loop, aC, DFG, FRET)

    .. todo :: Use kwargs with sensible defaults instead of relying only on positional arguments.


    """
    import mdtraj as md
    import numpy as np
    import pandas as pd

    topology = traj.topology
    coord = traj.xyz

    # JG debug
    #table, bonds = topology.to_dataframe()
    #atoms = table.values
    #print(atoms)

    # get the array of atom indices for the calculation of:
    #       * 7 ditances (a 7*2 array where each row contains indices of the two atoms for each distance)
    dis = np.zeros(shape=(2, 2), dtype=int, order='C')

    # name list of the features
    #feature_names = list()

    # parse the topology info
    '''
    The coordinates are located by row number (usually is atom index minus one, which is also why it's zero-based)
    by mdtraj but when the atom indices are not continuous there is a problem so a safer way to locate the coordinates
    is through row number (as a fake atom index) in case the atom indices are not continuous.
    '''

    ### distances
    if topology.select(f"chainid 0 and residue {key_res['group3'][3]} and name CA"):
        dis[0][0] = topology.select(f"chainid 0 and residue {key_res['group3'][3]} and name CA")
    if topology.select(f"chainid 0 and residue {key_res['group3'][2]} and name CZ"):
        dis[0][1] = topology.select(f"chainid 0 and residue {key_res['group3'][2]} and name CZ")
    if topology.select(f"chainid 0 and residue {key_res['group2'][3]} and name CA"):
        dis[1][0] = topology.select(f"chainid 0 and residue {key_res['group2'][3]} and name CA")
    dis[1][1] = dis[0][1]
    #feature_names.append('Dunbrack_D1')
    #feature_names.append('Dunbrack_D2')

    # check if there is any missing coordinates; if so, skip dihedral/distance calculation for those residues
    check_flag = 1

    for i in range(len(dis)):
        if 0 in dis[i]:
            dis[i] = [0,0]
            check_flag = 0
    #if check_flag:
        #print("There is no missing coordinates.  All dihedrals and distances will be computed.")


    distances = md.compute_distances(traj, dis)

    # clean up
    del traj, dis
    return distances
