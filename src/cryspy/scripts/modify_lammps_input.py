#!/usr/bin/env python3

### Header ###

# This script modifies the LAMMPS input file
# It requires the file calculated by LAMMPS-Interface

# It is included in stage 2

### Import modules ###

### Main ###

filename = "in.p1_structure"

in_data = open(filename, 'r').readlines()

read_data_index = 0
min_style_index = 0

for line in in_data:

    split_line = line.split()
    
    if len(split_line) > 0:
    
        if split_line[0] == "read_data":
        
            read_data_index = in_data.index(line)
            
for line in in_data:

    split_line = line.split()
    
    if len(split_line) > 0:
            
        if split_line[0] == "####":
        
            min_style_index = in_data.index(line)
                        
in_data[read_data_index] = "read_data   all_input_data\n"

line_2 = "dump_modify    structure_xyzmov element C H N O\n"
line_1 = "dump           structure_xyzmov 1 xyz 100000 structure_mov.xyz\n"

in_data.insert((min_style_index + 2), line_2)
in_data.insert((min_style_index + 2), line_1)

in_data.append("undump   structure_xyzmov\n\n")
in_data.append("print \" \"\n")
in_data.append("print \"----------------\"\n")
in_data.append("print \"Cell Data\"\n")
in_data.append("print \"----------------\"\n")
in_data.append("print \" \"\n\n")
in_data.append("variable x0 equal xlo\n")
in_data.append("variable x1 equal xhi\n")
in_data.append("variable y0 equal ylo\n")
in_data.append("variable y1 equal yhi\n")
in_data.append("variable z0 equal zlo\n")
in_data.append("variable z1 equal zhi\n")
in_data.append("variable XY equal xy\n")
in_data.append("variable XZ equal xz\n")
in_data.append("variable YZ equal yz\n\n")
in_data.append("print \"x0 = ${x0}\"\n")
in_data.append("print \"x1 = ${x1}\"\n")
in_data.append("print \"y0 = ${y0}\"\n")
in_data.append("print \"y1 = ${y1}\"\n")
in_data.append("print \"z0 = ${z0}\"\n")
in_data.append("print \"z1 = ${z1}\"\n")
in_data.append("print \"XY = ${XY}\"\n")
in_data.append("print \"XZ = ${XZ}\"\n")
in_data.append("print \"YZ = ${YZ}\"\n")

out_filename = "lammps.IN_1"

f_out = open(out_filename, "w")
f_out.writelines(in_data)
f_out.close()

### End ###
