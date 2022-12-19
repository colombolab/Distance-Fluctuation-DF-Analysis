# Distance-Flusctuation (DF) Analysis

What DF matrices are?
The analysis of the results of a MD simulation can performed using the Distance Fluctuation matrices (DF), based on the Coordination Propensity hypothesis:

![CodeCogsEqn-58-2](https://user-images.githubusercontent.com/70646674/208474307-0deb70ce-255e-47fc-8b92-cfecb1ab873d.png)


low CP values, corresponding to low pair-distance fluctuations, highlight groups of residues that move in a mechanically coordinated way.

# How to use the script

• Requisites

  - Python 3.0 (or newer version)

  - Numpy

  - Scipy

• Usage

The script can analyze a MD trajectory and identify the coordinated motions between residues. It can then filter the output matrix based on the distance to identify long-range coordinated motions.

The script can work both using only C-alphas (using either a pdb or a xyz file) or the sidechains (using a pdb file).

For more information run:

python3 distance_fluctuation.py -h

optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  -ext {pdb,xyz}, --file_ext {pdb,xyz}
                        Input trajectory file iformat (options: .pdb or .xyz)
  -n N_AA, --n_aa N_AA  Number of aminoacids in each frame of the .pdb
                        trajectory

Optional arguments:
  -i IN_NAME, --in_name IN_NAME
                        Name of the trajectory file (default: trj)
  -s S_FRAME, --s_frame S_FRAME
                        Index of the initial frame (default = 0)
  -e E_FRAME, --e_frame E_FRAME
                        Index of the last frame (default = end)
  -c CUTOFF, --cutoff CUTOFF
                        Cutoff for the distance fluctuation analysis (default = 5)
  -t TOLERANCE, --tolerance TOLERANCE
                        Tolerance for the distance fluctuation analysis (default = 0)
  -p {all,c-a}, --pdb_type {all,c-a}
                        Type of pdb file submitted (only C-alpha, all protein atoms)
  -l {s,seq,v,volume}, --local {s,seq,v,volume}
                        Type of local cutoff used, sequence or volume (default = seq)
  -r RADIUS, --radius RADIUS
                        Radius of the area cutoff (default = 7A)
  -rs RES_START, --res_start RES_START
                        Starting residue if the PDB does not start with res.num. 1
  -ro RMS_OUT, --rms_out RMS_OUT
                        RMS distance output filename (without extension)
  -ao AVG_OUT, --avg_out AVG_OUT
                        Average distance output filename (without extension)
  -so SEQ_OUT, --seq_out SEQ_OUT
                        Profile sequence output filename (without extension)
  -da DIST_AN, --dist_an DIST_AN
                        Distance analysis output filename (without extension)
  --blocks BLOCKS       File with domains borders

  -rb RMS_B_OUT, --rms_b_out RMS_B_OUT
                        RMS distance blocks output filename

# Read the output

The script generate different output file.

A) Average Distance (avgdist) file:

  The name of the file is avgdist_out_x_y_type.dat with x = start frame of the analysis, y = end frame of the analysis, type = type of analysis (CA or sidechains)

  The file contains a matrix using the residue indexes as axes and the average value of the distance between the residues as the data (r1 r2 avgdist).

  The distance is calculated as the average of the euclidean distance between the residues.

B) Distance Fluctuation (rmsdist) file:

  The name of the file is rmsdist_out_x_y_type.dat with x = start frame of the analysis, y = end frame of the analysis, type = type of analysis (CA or sidechains).

  The file contains a matrix using the residue indexes as axes and the distance fluctuation between the residues as the data (r1 r2 rmsdist).

  The distance fluctuation is calculated for the residues that are at least 3 residues away from each other (x-2 to x+2) as follow:

  1) Calculate the average euclidean distance between the residues (either CA or center of mass)

  2) Calculate the average distance vector

  3) Substract the distance fluctuation to the average distance

  4) Calculate the power of the difference between distance fluctuation and local fluctuation

  5) Filter the values of the close residues (1-x, x, 1+x)

  6) Divide the obtained value for the number of frames

C) Profile Sequence (profile_sequence) file:

  The name of the file is profile_sequence_x_y_type.dat with x = start frame of the analysis, y = end frame of the analysis, type = type of analysis (CA or sidechains).

  The file is a vector containing the residue number and the local flucutation value.

  The local fluctuation is calculated as the average fluctuation of the residues close to each other (that is the residues ranging from x-2 to x+2 with x = residue number).

  Since the value contains the average distance fluctuation for a range of residues the output starts from residue 4 and ends at residue n-3).

D) Distance Analysis (distance_analysis_out) file:

  The name of the file is distance_analysis_out_c_x_y_type.dat with c = cutoff value, x = start frame of the analysis, y = end frame of the analysis, type = type of analysis (CA or sidechains).

  This fail contains the distance fluctuation filtered by the cutoff value (the value are kept if the distance fluctuation is smaller than the cutoff value).

  The cutoff value is calculated as the sum between the nearcutoff and the tolerance value specified by the user.

  The nearcutoff can be calcolated either using sequence proximity or 3d proximity and can be specified by command line using the option -l or --local (see the command line help for further details).

  In the case of sequence proximity it is calculated as the average of the distance fluctuation value for residues in the range of x-2,x+2 (x = residue index)

  In the case of 3d proximity is calculated as the average distance fluctuation of the residues within a certain radius from the current one (default value = 6.5 A)

E) Blocks Averaging: (rmsdist_out_b)

It is possible to obtain a distance fluctuation matrix average on protein domains ( or blocks).

The script requires an input file with the borders of the protein blocks defined by the user, in the form of two columns: the first with the lower limit and the second with the upper limit.

# References
Morra, G.; Potestio, R.; Micheletti, C.; Colombo, G., Corresponding Functional Dynamics across the Hsp90 Chaperone Family: Insights from a Multiscale Analysis of MD Simulations, PLOS Computational Biology 8(3): e1002433. https://doi.org/10.1371/journal.pcbi.1002433
