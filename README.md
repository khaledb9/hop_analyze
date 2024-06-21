# hop_analyze

This code works side-by-side with https://github.com/danis-b/HOPPING

After extracting the hopping parameters from wannier90_hr.dat file: 

1- Use the label code with out.dat to include proper labelings and numbering for your material.

2- Execute the output unix commands to create .dat hopping files for each pair in your system.

3- Use the print_hop code to print useful information about each pair hopping. 


Note: the value of hopping in the columns corresponds to the highest value in the matrix 

Example Usage: 
python print_hop.py hop.Mo-Mo.dat 

Example Output:

Atoms Pair  |   Max Abs Value (Hopping)  |  Phase   |     Group     |  Radius
