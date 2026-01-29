#!/usr/bin/env python3

### Header ###

# This script reads the LAMMPS output from a CrySPY run
# It calculates densities from the molecular weight and the cell volume
# It generates density vs relative energy plots

### Import modules ###

from matplotlib import pyplot

import os

### Create tools ###

def isfloat(num):
    
    try:
        float(num)
        return True
    except ValueError:
        return False

### Read minimised energies and cell volumes ###

# Get user input filename
# This should be the logfile specified in LAMMPS interface

print("Please enter name of logfile: ")
logfile = input()

print("Please enter Molecular Weight of Unit Cell: ")
mr = float(input())

print("Please enter the number of atoms in the unit cell")
num_of_atoms = int(input())

print("Please enter the number types of atoms in the unit cell")
type_of_atoms = int(input())

print("Please enter number of structures generated: ")
struc_ctr = int(input())

data_rho = []
data_en = []

for entry in range(0, struc_ctr):

    directory = 'work/%i/'%entry
    
    ofile = directory + logfile

    print(ofile) # For debugging
    
    output = open(ofile, 'r').readlines()
    
    loop_ctr = len(output)

    vol = -1.0
    tot_en = -1.0

    trash = 1.0 # Don't ask
    
    for line in range(0, loop_ctr):
        # find the last density entry in the file

       if len(output[line].split()) != 7:

            trash = 1.1

       else:
	
            if output[line].split()[6] == 'Volume':
                vol_line = output[line + 2]

                if isfloat(vol_line.split()[6]) == True:

                    vol = float(vol_line.split()[6])
                    tot_en = float(vol_line.split()[4])
    
    # Calculate density and convert from Mr/A3 to g/cm3
    
    if vol != -1.0:         
        
        rho = (mr / vol) * 1.660578
    
        # Calculate energy per atom and convert from kcal/mol to eV
    
        atom_en = (tot_en / 23.06) / num_of_atoms

        data_rho.append(rho)
        data_en.append(atom_en)

### Plot graphs ###

ens = []

data_en_ctr = len(data_en)

print("Number of energies\n") # For debugging
print(data_en_ctr) # For debugging
print(data_en) # For debugging
print("\n") # For debugging0

for entry_en in range(0, data_en_ctr):

    ens.append(data_en[entry_en])

# Get user input for Crystal Structure bechmark

print("Do you want absolute energy (yes/no)?: ")

abs_en_request = input()

if abs_en_request == "no":

    print("Please provide reference energy: ")

    minen = float(input())

    ens = [en_save - minen for en_save in ens]
    
else:

    minen = min(ens)
    
    ens = [en_save - minen for en_save in ens]

print("The ens array") # For debugging
print(ens) # For debugging

# This for loop is designed to remove outliers
# Modify value as applicable
# This section has been rewritten due to a floating point precision error

rem_ctr = len(ens)

trash = []

for result in range(0, rem_ctr):

    if ens[result] > 1.0:

        # rem = data_en.index(result)
        # del ens[rem]
        # del data_rho[rem]
        
        app_trash = result
        
        trash.append(app_trash)
        
print(trash)

trash_flag = 0

trash_ctr = len(trash)

for result_2 in range(0, trash_ctr):

	removal = trash[result_2] - trash_flag
        
	del ens[removal]
	del data_rho[removal]
	
	trash_flag += 1 

print("Do you want absolute rho (yes/no)?: ")

abs_request = input()

data_rho_ctr = len(data_rho)

print("Number of densities\n") # For debugging
print(data_rho_ctr) # For debugging
print(data_rho) # For debugging
print("\n") # For debugging

if abs_request == "yes":

   rhos = []
   
   for entry_rho in range(0, data_rho_ctr):
   
       rhos.append(data_rho[entry_rho])
	
else:

   rhos = []
   
   for entry_rho in range(0, data_rho_ctr):
   
       rhos.append(data_rho[entry_rho])
   
   print("Please provide reference density: ")
   
   minrho = float(input())
   
   rhos = [rho_save - minrho for rho_save in rhos]

# spg = [2 * entry['Sgp_num'] for entry in data]

print("The rhos array") # For debugging
print(rhos) # For debugging

# This for loop is designed to remove outliers
# Modify value as applicable

#for result in ens:
#
#   if result > 1.0:
#
#       rem = ens.index(result)
#
#      del ens[rem]
#      del rhos[rem]
 
pyplot.scatter(ens, rhos, c="Blue", alpha=0.5)
pyplot.xlabel('Rel. energy per atom (eV)')
   
if abs_request == "yes":   
   
   pyplot.ylabel('Density (g / cm3)')
   
else:

   pyplot.ylabel('Rel. density (g / cm3)')

#pyplot.show()
pyplot.savefig('rho_en.png')

print("Density vs Energy plot generated.\n")


print("Would you like to perform single point DFT for all structures? (yes/no)")

create_sp_dft = input()

if create_sp_dft == "yes":

    print("Please enter the name of the structural output from LAMMPS")
    struc_out = input()

    print("Please input the desired k point grid(# # #):")
    k_points = input()

    for low_dir in range(0, struc_ctr):
        print(low_dir)
        # Work in this loop so there's only one structure file open at a time
        
        # Get the target results directory

        target_dir = 'work/%i/'%low_dir
        
        # Create the DFT directory
        
        new_dir = target_dir + "DFT"
        
        os.mkdir(new_dir)
        
        # Get the structure output from the LAMMPS calculation
        
        struc_file = target_dir + struc_out
        
        ### Write the Quantum Espresso input ###
        
        # Use readlines to get box
        
        log_dir = target_dir + logfile
        
        cell_out = open(log_dir, 'r').readlines()
        log_ctr = len(cell_out)
        
        xlo = -1.0
        xhi = -1.0
        ylo = -1.0
        yhi = -1.0
        zlo = -1.0
        zhi = -1.0
        
        xy = -1.0
        xz = -1.0
        yz = -1.0
        try:
            xlo = float(cell_out[-18].split()[2])
            xhi = float(cell_out[-16].split()[2])
            ylo = float(cell_out[-14].split()[2])
            yhi = float(cell_out[-12].split()[2])
            zlo = float(cell_out[-10].split()[2])
            zhi = float(cell_out[-8].split()[2])
        
            xy = float(cell_out[-6].split()[2])
            xz = float(cell_out[-4].split()[2])
            yz = float(cell_out[-2].split()[2])
        except:
            print("Job didn't finish")
            continue

        if xlo == -1.0:
        
            print("WARNING: xlo not found")
            
        if xhi == -1.0:
        
            print("WARNING: xhi not found")
            
        if ylo == -1.0:  
        
            print("WARNING: ylo not found")
            
        if yhi == -1.0:
        
            print("WARNING: yhi not found")
            
        if zlo == -1.0:
        
            print("WARNING: zlo not found")
            
        if zhi == -1.0:
        
            print("WARNING: zhi not found")
            
        if xy == -1.0:
        
            print("WARNING: xy not found")
            
        if xz == -1.0:
        
            print("WARNING: xz not found")
            
        if yz == -1.0:
        
            print("WARNING: yz not found")
        
        # Calculate box vectors
        
        a1 = xhi - xlo
        a2 = 0.0
        a3 = 0.0
        
        b1 = xy
        b2 = yhi - ylo
        b3 = 0.0
        
        c1 = xz
        c2 = yz
        c3 = zhi - zlo
        
        vec_1 = str(a1) + " " + str(a2) + " " + str(a3)
        vec_2 = str(b1) + " " + str(b2) + " " + str(b3)
        vec_3 = str(c1) + " " + str(c2) + " " + str(c3)
        
        # Use readlines to get coordinates
        
        coords_out = open(struc_file, 'r').readlines()
        
        coordinates = []
        
        num_of_atoms = int(num_of_atoms)
        
        for coord_atom in range(0, num_of_atoms):
        
            tkr = -1 - coord_atom
            
            coordinates.append(coords_out[tkr])
            
        num_of_coords = len(coordinates)
        
        if num_of_coords != num_of_atoms:
        
            print("WARNING! Error in atom coordinate import!")
        
        # Create QE input file
        
        input_file = new_dir + "/pwscf.in"
        
        with open(input_file, 'a+') as in_f:
        
            # Write in the input file options
             # Manually modify don't forget!!

            in_f.write("&control\n")
            
            in_f.write("   title = \'CHNO\'\n")
            in_f.write("   calculation = \'scf\'\n")
            in_f.write("   nstep = 1000\n")
            in_f.write("   restart_mode = \'from_scratch\'\n")
            in_f.write("   pseudo_dir = \'/home/cm/cmpm4/Pseudo/SSSP_1_3_0_PBE_efficiency\'\n")
            in_f.write("   outdir = \'./outdir/\'\n")
            in_f.write("   etot_conv_thr = 1.0D-5\n")
            in_f.write("   forc_conv_thr = 1.0D-5\n")           
            

            in_f.write("/\n\n") # Include blank line
            
            in_f.write("&system\n")
            
            in_f.write("   ibrav = 0\n")
            in_f.write(f"   nat = {num_of_atoms}\n") # Modify number of atoms here
            in_f.write(f"   ntyp = {type_of_atoms}\n") # Modify number of atom types here
            in_f.write("   ecutwfc = 60.00\n")
            in_f.write("   ecutrho = 480.00\n")
            in_f.write("   occupations = \'smearing\'\n")
            in_f.write("   degauss = 0.01\n")
            in_f.write("   vdw_corr = \'grimme-d3\'\n")
            in_f.write("   dftd3_version = 6\n")
            
            in_f.write("/\n\n") # Include blank line
            
            in_f.write("&electrons\n")
            in_f.write("conv_thr = 1.D-6\n")
            in_f.write("mixing_beta = 0.5D0\n")
            in_f.write("diagonalization = \'david\'\n")
            in_f.write("/\n\n") # Include blank line
            
            
            # Atomic Species will need to be changed if the atom species change
            in_f.write("ATOMIC_SPECIES\n")
            
            in_f.write("   H  1.008  H.pbe-rrkjus_psl.1.0.0.UPF\n")
            in_f.write("   C  12.011 C.pbe-n-kjpaw_psl.1.0.0.UPF\n")
            in_f.write("   N  14.00  N.pbe-n-radius_5.UPF\n")
            in_f.write("   O  16.00  O.pbe-n-kjpaw_psl.0.1.UPF\n")
            in_f.write("\n") # Blank line included here to make atom species change easier
            
            # Print the lattice parameters here
            
            in_f.write("CELL_PARAMETERS angstrom\n")
            in_f.write(vec_1)
            in_f.write("\n")
            in_f.write(vec_2)
            in_f.write("\n")
            in_f.write(vec_3)
            in_f.write("\n")
            
            # Print the atomic positions here
            
            in_f.write("ATOMIC_POSITIONS angstrom\n")
            
            for print_ctr in range(0, num_of_atoms):
            
                in_f.write(coordinates[print_ctr])
                # in_f.write("\n")
            
            # Don't forget the k-points, put them here
            # Have user input for these
            
            # print("Please input the desired k point grid(# # #):")
            # k_points = input()
            in_f.write("\nK_POINTS automatic\n")
            kp_input = k_points + " 0 0 0"
            in_f.write(kp_input)
            
            in_f.close()
            
        # Write BASH file
        

print("Thank you for choosing us to generate your output, have a nice day!")
print("Would you like to perform vc-relax DFT minimisation? (yes/no)")

create_dft = input()

if create_dft == "yes":

    print("How many of the top structures do you wish to minimise?")
    dft_ctr = int(input())

    # Apply all steps to both energy and density lists
    # Create new lists
    
    sorted_energies = ens.copy()
    sorted_densities = rhos.copy()
    
    # Sort new lists
    # Get lowest energies and largest densities
    
    sorted_energies.sort()
    sorted_densities.sort(reverse=True)
    
    # Search original lists for the indexes of the first 10 values
    
    low_ens = []
    low_rhos = []
    
    for low_en in range(0, dft_ctr):
    
        low_ens.append(sorted_energies[low_en])
    
    for low_rho in range(0, dft_ctr):
    
        low_rhos.append(sorted_densities[low_rho])
        
    # Use sorted lists to acquire indexes of entries in ens and rhos
    
    target_list = []
        
    for target_en in range(0, dft_ctr):
    
        target_en_index = ens.index(low_ens[target_en])
        
        target_flag = 0
        
        for target in target_list:
        
            if target_en_index == target:
            
                target_flag = target_flag + 1
                
        if target_flag == 0:
        
            target_list.append(target_en_index)
            
    for target_rho in range(0, dft_ctr):
    
        target_rho_index = rhos.index(low_rhos[target_rho])
        
        target_flag_rho = 0
        
        for target_rho in target_list:
        
            if target_rho_index == target_rho:
            
                target_flag_rho = target_flag_rho + 1
                
        if target_flag_rho == 0:
        
            target_list.append(target_rho_index)
            
    # Create list of directories from the index list
    
    dir_ctr = len(target_list)
    print(target_list)
    
    print("Please enter the name of the structural output from LAMMPS")
    struc_out = input()
    
    print("Please input the desired k point grid(# # #):")
    k_points = input()
    
    # Create path for BASH submission
    
    bash_file = "manual_dft_submission"
    
    # Loop connecting LAMMPS and Quantum Espresso
    
    for low_dir in range(0, dir_ctr):
        print(low_dir)
        # Work in this loop so there's only one structure file open at a time
        
        # Get the target results directory

        target_dir = 'work/%i/'%target_list[low_dir]
        
        # Create the DFT directory
        
        new_dir = target_dir + "DFT"
        
        os.mkdir(new_dir)
        
        # Get the structure output from the LAMMPS calculation
        
        struc_file = target_dir + struc_out
        
        ### Write the Quantum Espresso input ###
        
        # Use readlines to get box
        
        log_dir = target_dir + logfile
        
        cell_out = open(log_dir, 'r').readlines()
        log_ctr = len(cell_out)
        
        xlo = -1.0
        xhi = -1.0
        ylo = -1.0
        yhi = -1.0
        zlo = -1.0
        zhi = -1.0
        
        xy = -1.0
        xz = -1.0
        yz = -1.0
        try:
            xlo = float(cell_out[-18].split()[2])
            xhi = float(cell_out[-16].split()[2])
            ylo = float(cell_out[-14].split()[2])
            yhi = float(cell_out[-12].split()[2])
            zlo = float(cell_out[-10].split()[2])
            zhi = float(cell_out[-8].split()[2])
        
            xy = float(cell_out[-6].split()[2])
            xz = float(cell_out[-4].split()[2])
            yz = float(cell_out[-2].split()[2])
        except:
            print("Job didn't finish")
            continue



#        for log_line in range(0, log_ctr):
#        
#            if bool(cell_out[log_line]):
#            
#                if cell_out[log_line].split()[0] == "x0":
#            
#                    xlo = cell_out[log_line].split()[2]
#                
#                elif cell_out[log_line].split()[0] == "x1":
#                
#                    xhi = cell_out[log_line].split()[2]
#                    
#                elif cell_out[log_line].split()[0] == "y0":
#                
#                    ylo = cell_out[log_line].split()[2]
#                    
#                elif cell_out[log_line].split()[0] == "y1":
#                
#                    yhi = float(cell_out[log_line].split()[2])
#                    
#                elif cell_out[log_line].split()[0] == "z0":
#                
#                    zlo = float(cell_out[log_line].split()[2])
#                    
#                elif cell_out[log_line].split()[0] == "z1":
#                
#                    zhi = float(cell_out[log_line].split()[2])
#                    
#                elif cell_out[log_line].split()[0] == "XY":
#                
#                    xy = float(cell_out[log_line].split()[2])
#                    
#                elif cell_out[log_line].split()[0] == "XZ":
#                
#                    xz = float(cell_out[log_line].split()[2])
#                    
#                elif cell_out[log_line].split()[0] == "YZ":
#                
#                    yz = float(cell_out[log_line].split()[2])
            
        if xlo == -1.0:
        
            print("WARNING: xlo not found")
            
        if xhi == -1.0:
        
            print("WARNING: xhi not found")
            
        if ylo == -1.0:  
        
            print("WARNING: ylo not found")
            
        if yhi == -1.0:
        
            print("WARNING: yhi not found")
            
        if zlo == -1.0:
        
            print("WARNING: zlo not found")
            
        if zhi == -1.0:
        
            print("WARNING: zhi not found")
            
        if xy == -1.0:
        
            print("WARNING: xy not found")
            
        if xz == -1.0:
        
            print("WARNING: xz not found")
            
        if yz == -1.0:
        
            print("WARNING: yz not found")
        
        # Calculate box vectors
        
        a1 = xhi - xlo
        a2 = 0.0
        a3 = 0.0
        
        b1 = xy
        b2 = yhi - ylo
        b3 = 0.0
        
        c1 = xz
        c2 = yz
        c3 = zhi - zlo
        
        vec_1 = str(a1) + " " + str(a2) + " " + str(a3)
        vec_2 = str(b1) + " " + str(b2) + " " + str(b3)
        vec_3 = str(c1) + " " + str(c2) + " " + str(c3)
        
        # Use readlines to get coordinates
        
        coords_out = open(struc_file, 'r').readlines()
        
        coordinates = []
        
        num_of_atoms = int(num_of_atoms)
        
        for coord_atom in range(0, num_of_atoms):
        
            tkr = -1 - coord_atom
            
            coordinates.append(coords_out[tkr])
            
        num_of_coords = len(coordinates)
        
        if num_of_coords != num_of_atoms:
        
            print("WARNING! Error in atom coordinate import!")
        
        # Create QE input file
        
        input_file = new_dir + "/pwscf.in"
        
        with open(input_file, 'a+') as in_f:
        
            # Write in the input file options
             # Manually modify don't forget!!

            in_f.write("&control\n")
            
            in_f.write("   title = \'CHNO\'\n")
            in_f.write("   calculation = \'vc-relax\'\n")
            in_f.write("   nstep = 1000\n")
            in_f.write("   restart_mode = \'from_scratch\'\n")
            in_f.write("   pseudo_dir = \'/home/c/cmpm4/Pseudo/SSSP\'\n")
            in_f.write("   outdir = \'./outdir/\'\n")
            in_f.write("   etot_conv_thr = 1.0D-5\n")
            in_f.write("   forc_conv_thr = 1.0D-5\n")           
            

            in_f.write("/\n\n") # Include blank line
            
            in_f.write("&system\n")
            
            in_f.write("   ibrav = 0\n")
            in_f.write(f"   nat = {num_of_atoms}\n") # Modify number of atoms here
            in_f.write(f"   ntyp = {type_of_atoms}\n") # Modify number of atom types here
            in_f.write("   ecutwfc = 60.00\n")
            in_f.write("   ecutrho = 480.00\n")
            in_f.write("   occupations = \'smearing\'\n")
            in_f.write("   degauss = 0.01\n")
            in_f.write("   vdw_corr = \'grimme-d3\'\n")
            in_f.write("   dftd3_version = 6\n")
            
            in_f.write("/\n\n") # Include blank line
            
            in_f.write("&electrons\n")
            in_f.write("conv_thr = 1.D-6\n")
            in_f.write("mixing_beta = 0.5D0\n")
            in_f.write("diagonalization = \'david\'\n")
            in_f.write("/\n\n") # Include blank line
            
            in_f.write("&ions\n")
            
            in_f.write("ion_dynamics = \'bfgs\'\n")
            
            in_f.write("/\n\n") # Include blank line
            
            in_f.write("&cell\n")
            
            in_f.write("cell_dynamics = \'bfgs\'\n")
            in_f.write("press_conv_thr = 0.1\n")
            in_f.write("cell_factor = 3.0\n")
            
            in_f.write("/\n\n") # Include blank line
            
            # Atomic Species will need to be changed if the atom species change
            in_f.write("ATOMIC_SPECIES\n")
            
            in_f.write("   H  1.008  H.pbe-rrkjus_psl.1.0.0.UPF\n")
            in_f.write("   C  12.011 C.pbe-n-kjpaw_psl.1.0.0.UPF\n")
            in_f.write("   N  14.00  N.pbe-n-radius_5.UPF\n")
            in_f.write("   O  16.00  O.pbe-n-kjpaw_psl.0.1.UPF\n")
            in_f.write("\n") # Blank line included here to make atom species change easier
            
            # Print the lattice parameters here
            
            in_f.write("CELL_PARAMETERS angstrom\n")
            in_f.write(vec_1)
            in_f.write("\n")
            in_f.write(vec_2)
            in_f.write("\n")
            in_f.write(vec_3)
            in_f.write("\n")
            
            # Print the atomic positions here
            
            in_f.write("ATOMIC_POSITIONS angstrom\n")
            
            for print_ctr in range(0, num_of_atoms):
            
                in_f.write(coordinates[print_ctr])
                # in_f.write("\n")
            
            # Don't forget the k-points, put them here
            # Have user input for these
            
            # print("Please input the desired k point grid(# # #):")
            # k_points = input()
            in_f.write("\nK_POINTS automatic\n")
            kp_input = k_points + " 0 0 0"
            in_f.write(kp_input)
            
            in_f.close()
            
        # Write BASH file
        
        with open(bash_file, 'a+') as bash_f:
        
            change_dir = "cd " + new_dir
            bash_f.write(change_dir)
            bash_f.write("\n")
            cp_cmd = "cp /mnt/gpfs01/home/cm/cmpm4/git/CrySPY/qe_run ."
            bash_f.write(cp_cmd)
            bash_f.write("\n")
            bash_f.write("qsub qe_runscript_10\n")
            bash_f.write("cd ../../../../\n\n")
            
            bash_f.close()

print("Thank you for choosing us to generate your output, have a nice day!")

### End ###
