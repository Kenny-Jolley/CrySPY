#!/usr/bin/env python3

### Header ###

# This script writes a lammps runscript
# !!! THIS IS FOR LOVELACE !!!

# It is included in stage 2

### Import modules ###

### Main ###

run_text = []

run_text.append("#!/bin/bash\n")
run_text.append("#PBS -N LAMMPS-eam-N1\n")
run_text.append("#PBS -l nodes=1:ppn=1\n")
run_text.append("#PBS -l walltime=100:00:00\n")
run_text.append("#PBS -A Jolley2022a\n")
run_text.append("#PBS -q compute\n")
run_text.append("   \n")
run_text.append("module purge\n")
run_text.append("module load LAMMPS/23Jun2022-intel-2022b\n")
run_text.append("   \n")
run_text.append("export NCPUS=`wc -l < $PBS_NODEFILE`\n")
run_text.append("export OMP_NUM_THREADS=2\n")
run_text.append("   \n")
run_text.append("cd $PBS_O_WORKDIR\n")
run_text.append("   \n")
run_text.append("modify_potential.py\n")
run_text.append("prepare_cat.py\n")
run_text.append("cat lammps.lattice.dat data.potentials > all_input_data\n")
run_text.append("   \n")
run_text.append("mpirun -np $NCPUS lmp_mpi -sf intel -pk intel 0 omp 2 -in lammps.IN -log lammps.OUT\n")
run_text.append("   \n")
run_text.append("if [ -e \"CRASH\" ] ; then\n")
run_text.append("      sed -i -e \'3s/^.*$/skip/\' stat_job\n")
run_text.append("      exit 1\n")
run_text.append("fi\n")
run_text.append("   \n")
run_text.append("sed -i -e \'3s/^.*$/done/\' stat_job\n")
run_text.append("   \n")

# Write and save the text to file

run_file = "lammps_runscript"
    
f_run = open(run_file, "w")
f_run.writelines(run_text)
f_run.close()

### END ###
