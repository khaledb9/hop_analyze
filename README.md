# hop_analyze

This code works side-by-side with https://github.com/danis-b/HOPPING
After extracting the hopping parameters from wannier90_hr.dat file: 
1- Use the label code with out.dat to include proper labelings for your material and numbering.
2- Execute the output unix commands to create .dat hopping files for each pair in your system .
3- Use the print_hop code to print useful information about each pair hopping. 

Note: the value of hopping in the columns corresponds to the highest value in the matrix 

Example Usage: 
python print_hop.py hop.Mo-Mo.dat 

Example Output:
Atoms Pair     Max Abs Value (Hopping)    Phase              Group       Radius
Mo1-Mo2         0.0005                    -0.7952π             1NN    3.207019
Mo2-Mo1         0.0005                    -0.8789π             1NN    3.207019
Mo3-Mo5         0.5134                    -0.0075π             1NN    3.201611
Mo3-Mo4         0.2646                    -0.9976π             1NN    3.257957
Mo4-Mo5         0.4961                    -0.9982π             1NN    3.224519
Mo4-Mo3         0.3806                    -0.0066π             1NN    3.257957
Mo5-Mo3         0.4629                    -0.0199π             1NN    3.201611
Mo5-Mo4         0.4043                    -0.9950π             1NN    3.224519
Mo5-Mo8         0.4340                    -0.0045π             1NN    3.298487
Mo6-Mo7         0.0009                    -0.8858π             1NN    3.209084
Mo7-Mo6         0.0011                    -0.6872π             1NN    3.209084
Mo8-Mo9         0.0003                    -0.8976π             1NN    3.256762
Mo8-Mo5         0.4160                    -0.0030π             1NN    3.298487
Mo9-Mo8         0.0008                    -0.1652π             1NN    3.256762
Mo2-Mo9         0.0003                    -0.8976π             2NN    5.695201
Mo2-Mo5         0.0575                    -0.9505π             2NN    5.732545
Mo4-Mo7         0.0001                    -1.0000π             2NN    5.626167
