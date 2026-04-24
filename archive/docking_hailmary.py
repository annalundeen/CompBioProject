

# I have muscle proteins ABLIM1.pdb, MBYC3.pdb, and MYL2.pdb and toxin LqhIII.pdb in the folder CompBioProject/StructureVisualization

#there are several different motifs I want to test for each muscle protein

# want to make this whole pipeline test each unique protein + motif + toxin combination

#also need to figure out virtual environment to conda environment transition 



protein_motifs = {'MYL2': ['IDEMIK', 'TILN', 'IDEMIKE', 'NEEIDEMIKEAPG'], 'MYBPC': ['SFVPEGFAC', 'GQEIQMSGSKYIF', 'GQALAELIVQEKK', 'NFDLIQ', 'LAELIV'], 'ABLIM1': ['KVCGCDLAQGGFF']}


from Bio.PDB import PDBParser, PDBIO

parser = PDBParser()
structure = parser.get_structure("protein", "muscle_protein.pdb")

#protein
io = PDBIO()
io.set_structure(structure)
io.save("muscle_clean.pdb")

#toxin
structure = parser.get_structure("toxin", "toxin.pdb")
io.set_structure(structure)
io.save("toxin_clean.pdb")

###############
# command line
################
# need to install autodock tools
conda install -c conda-forge autodocktools

prepare_receptor4.py -r muscle_clean.pdb -o muscle.pdbqt
prepare_ligand4.py -l toxin_clean.pdb -o toxin.pdbqt

##############
# PYMOL
###############

#get motif binding center in pymol
select motif, pepseq KVCGCDLAQGGFF
center motif
get_position # [  22.426,  10.569, -41.382]
#example output: (12.3, -5.8, 22.1) #put this in config file

#######
# insert python script to create vina config files for each unique muscle protein + toxin + motif combination

#ex: ABLIM1, LqhIII, KVCGCDLAQGGFF, positon [22.426,  10.569, -41.382]
#would need to extract center x (first item in postion), center y (second item in positon), center z (third item in position)

####################
# .TXT file? Get a script to make these files?
####################
#create vina config file- example (do for each combo and get center)
receptor = muscle.pdbqt
ligand = toxin.pdbqt

center_x = 12.3
center_y = -5.8
center_z = 22.1

size_x = 20
size_y = 20
size_z = 20

exhaustiveness = 8 #increase to 32 for protein-protein interaction?
num_modes = 10

####################
# TERMINAL
###################

#run docking
vina --config config.txt --out docked.pdbqt --log docking.log

#install openbable
conda install -c conda-forge openbabel

#convert output to pdb
obabel docked.pdbqt -O docked.pdb



####################
# PYMOL
###################
#final visualization in pymol
# pymol muscle_clean.pdb docked.pdb

# show cartoon
# show sticks, docked
# color red, docked

#to highlight motif
# select motif, resi 45+46+47+48
# color yellow, motif