import mdtraj as md
import sys
from simtk.openmm.app import PDBFile
from features import featurize
from features import dunbrack_cluster

'''
# negtive example
pdbid = '4YNE'
chain = 'A'
flip = 0 
top_file = '/home/guoj1/data_projects/NTRK/NTRK1/apo_dfg/01.min_equi/receptor_equilibrated.pdb'
traj_file = '/home/guoj1/data_projects/NTRK/NTRK1/apo_dfg/02.sams_dunbrack/NTRK1_apo_4fs_traj.dcd'
'''
# test example (state 7)
pdbid = '1M17'
chain = 'A'
traj = md.load_pdb(f'{pdbid}_chain{chain}.pdb')

# JG debug
#traj = md.load('test_systems/2JIU_chainB.pdb') has been shown to be 0.8912438, 0.43516478
#table, bonds = traj.topology.to_dataframe()
#atoms = table.values
#print(atoms)

dihedrals, distances = featurize(pdb=f'{pdbid}', chain=f'{chain}', coord=traj, feature='conf')
print(distances)
# assign DFG states
assignment = dunbrack_cluster.assign(dihedrals, distances)

print(assignment)
