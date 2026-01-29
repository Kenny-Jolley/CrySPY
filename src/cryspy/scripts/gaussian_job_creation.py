#!/usr/bin/env python3

### Header ###

# This script creates Gaussian jobs to calculate partial charges
# It is included in stage 1

### Import modules ###

import os

### Main ###

# Get list of xyz files  

filename = "molelist"

moledata = open(filename, 'r').readlines()

molenames = []

for line in moledata:

    moledata_line = line.split()
    
    if len(moledata_line) > 0:
    
        molenames.append(moledata_line[1])
        
# molenames should be a list of molecule xyz names

# create Gaussian directories

dir_names = []

for molecule in molenames:

    molecule.rstrip()
    dir_name = molecule.rstrip(".xyz")
    
    os.mkdir(dir_name)
    
    dir_names.append(dir_name)
    
# create Gaussian input

for molecule in molenames:

    # extract the molecular coordinates from the xyz file

    coords = open(molecule, 'r').readlines()
    
    # modify the coordinates in the list to create a Gaussian input
    
    coords[0] = "chk=charges \n"
    
    coords[1] = "# PBEPBE/6-311G opt pop=full \n\n"
    coords.insert(2, "Title \n\n")
    
    # charge and multiplicity
    
    coords.insert(3, "0 1\n")

    # don't forget the blank line at the end of the file

    coords.append("\n")
    
    # create the input file and save the input
    
    tkr = molenames.index(molecule)
    
    f_name = dir_names[tkr] + "/charges.com"
    
    f = open(f_name, "w")
    f.writelines(coords)
    f.close()

# create a list for the runscript 
# !!! FOR LOVELACE ONLY !!!
# modify this section to change runscript input
# or modify manually in the file

run_text = []

run_text.append("#!/bin/bash\n")
run_text.append("#PBS -N g09\n")
run_text.append("#PBS -l nodes=1:ppn=20\n")
run_text.append("#PBS -l walltime=100:00:00\n")
run_text.append("#PBS -A Jolley2022a\n")
run_text.append("#PBS -q compute\n")
run_text.append("# vis, gpu\n")
run_text.append("   \n")
run_text.append("# Use less than 180GB mem\n")
run_text.append("   \n")
run_text.append("hostname\n")
run_text.append("   \n")
run_text.append("module load g09\n")
run_text.append("   \n")
run_text.append("cd $PBS_O_WORKDIR\n")
run_text.append("g09 charges.com charges.log\n")
run_text.append("   ")

for dir_path in dir_names:

    run_file = dir_path + "/PBS_g09"
    
    f_run = open(run_file, "w")
    f_run.writelines(run_text)
    f_run.close()
    
print("Gaussian file generation completed.")

### END ###
