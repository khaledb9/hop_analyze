# hop_analyze

This code works side-by-side with https://github.com/danis-b/HOPPING

After extracting the hopping parameters from wannier90_hr.dat file: 

1- Use the label code with out.dat from HOPPING.py to include proper labelings and numbering for your material. Adjust the label.py script for the number and type of elements.

2- Use the print_hop.py code to print the max abs value of each hopping block and filter the distances using the --gap options. 

3- A new feature is added to extract the on-site energies and their orbital characters in the script cf.py. You might want to use that only with atomic-like wannier functions (num_iter=0) for easier interpretation as the code does not modify the basis.
