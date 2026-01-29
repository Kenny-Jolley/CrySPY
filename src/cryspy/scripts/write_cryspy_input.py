#!/usr/bin/env python3

### Header ###

# This script generates the final cryspy.in file
# This should be the only one the user sees

# It is included in stage 2

### Import modules ###

### Main ###

# List of xyz file paths needed for input
# Number of each molecule is also needed

filename = "molelist"

moledata = open(filename, 'r').readlines()

molenames = []
molenums = []

for line in moledata:

    moledata_line = line.split()
    
    if len(moledata_line) > 0:
    
        molenames.append(moledata_line[1])
        molenums.append(moledata_line[0])

# Calculate number and type of atoms

atom_type_list = []

for molecule in molenames:

    coordfile = molecule
    coords = open(coordfile, 'r').readlines()
    
    col_num_var = int(molenames.index(molecule))
    num_mols = int(molenums[col_num_var])
    
    for mol_repeat in range(0, num_mols):
    
        for atom in coords:
        
            atom_line = atom.split()     
        
            if len(atom_line) > 1:
            
                atom_type_list.append(atom_line[0])
            
sort_atoms = sorted(atom_type_list)

num_of_atoms = len(sort_atoms)

atom_types = []
num_atoms_of_type = []

atom_types.append(sort_atoms[0])
type_ctr = 0
num_types = 1
comp_var = sort_atoms[0]

for atom in sort_atoms:

    if atom == comp_var:
        
        type_ctr = type_ctr + 1
        
    else:
    
        atom_types.append(atom)
        num_atoms_of_type.append(type_ctr)
        type_ctr = 1
        comp_var = atom

num_atoms_of_type.append(type_ctr)

# create CrySPY input

cryspy_text = []

cryspy_text.append("[basic]\n")
cryspy_text.append("algo = RS\n")
cryspy_text.append("calc_code = LAMMPS\n")
cryspy_text.append("tot_struc = 100\n")
cryspy_text.append("nstage = 1\n")
cryspy_text.append("njob = 250\n")
cryspy_text.append("jobcmd = qsub\n")
cryspy_text.append("jobfile = lammps_runscript\n")
cryspy_text.append("   \n")
cryspy_text.append("[structure]\n")
cryspy_text.append("struc_mode = mol\n")

# Total number of atoms

total_line = "natot = " + str(num_of_atoms) + "\n"

cryspy_text.append(total_line)

# Element list

ele_line = "atype = "

for a_type in atom_types:

    ele_line = ele_line + a_type + " "
    
ele_line = ele_line + "\n"

cryspy_text.append(ele_line)

# Number of each element

num_ele = "nat = "

for number in num_atoms_of_type:

    num_ele = num_ele + str(number) + " "

num_ele = num_ele + "\n"

cryspy_text.append(num_ele)

# XYZ file location

xyz_line = "mol_file = "

for molpath in molenames:

    molpath = "./" + molpath + "  "
    
    xyz_line = xyz_line + molpath

xyz_line = xyz_line + "\n"

cryspy_text.append(xyz_line)

# Number of each molecule

molnum_line = "nmol = "

for molnumb in molenums:

    molnumb = molnumb + " "
    
    molnum_line = molnum_line + molnumb
    
molnum_line = molnum_line + "\n"

cryspy_text.append(molnum_line)

# This should be modified by the user

cryspy_text.append("spgnum = 1\n")

# Minimum distance matrix

matrix_ctr = len(atom_types)

dist_vals = [1.0] * matrix_ctr

for x in range(0, matrix_ctr):

    # create string
    
    mindist_str = "mindist_" + str(x) + " = "
    
    # loop to add numerical values to string
    
    for y in range(0, matrix_ctr):
    
        mindist_str = mindist_str + str(dist_vals[y]) + " "
    
    mindist_str = mindist_str + "\n"
    
    # append string to cryspy_text
    
    cryspy_text.append(mindist_str)

# Write the rest of the CrySPY input

cryspy_text.append("symprec = 0.05\n")
cryspy_text.append("   \n")
cryspy_text.append("[LAMMPS]\n")
cryspy_text.append("lammps_infile = lammps.IN\n")
cryspy_text.append("lammps_outfile = lammps.OUT\n")
cryspy_text.append("lammps_potential = data.potentials\n")
cryspy_text.append("lammps_data = lammps.lattice.dat\n")
cryspy_text.append("   \n")
cryspy_text.append("[option]\n")

# save CrySPY input to file

cryspy_file = "cryspy.in"

f_out = open(cryspy_file, "w")
f_out.writelines(cryspy_text)
f_out.close()

print("CrySPY input completed")

### END ###
