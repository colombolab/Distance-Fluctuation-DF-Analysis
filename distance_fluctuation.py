#!/usr/bin/python

import argparse, re, sys, os, textwrap
from itertools import islice, zip_longest
import tempfile
import numpy as np
import math
from scipy.spatial import distance

masses_pdb = { 'C'  :  0.0000, 'CA' : 12.0107, 'CB' : 12.0107, 'CD' : 12.0107, 'CD1': 12.0107, 'CD2': 12.0107, 'CE' : 12.0107, 'CE1': 12.0107, 'CE2': 12.0107,
               'CE3': 12.0107, 'CG' : 12.0107, 'CG1': 12.0107, 'CG2': 12.0107, 'CH2': 12.0107, 'CZ' : 12.0107, 'CZ2': 12.0107, 'CZ3': 12.0107,

               'N'  :  0.0000, 'ND1': 14.0067, 'ND2': 14.0067, 'NE' : 14.0067, 'NE1': 14.0067, 'NE2': 14.0067, 'NH1': 14.0067, 'NH2': 14.0067, 'NZ' : 14.0067,

               'O'  :  0.0000, 'OD' : 15.9994, 'OD1': 15.9994, 'OD2': 15.9994, 'OE1': 15.9994, 'OE2': 15.9994, 'OG' : 15.9994, 'OG1': 15.9994, 'OH' : 15.9994, 'OXT': 0.0000,

               'SD' : 32.0650, 'SG' : 32.0650,

               'H'   : 0.00000, 'H1'  : 0.00000, 'H2'  : 0.00000, 'H3'  : 0.00000, 'HA'  : 1.00794, 'HA2' : 1.00794, 'HA3' : 1.00794, 'HB'  : 1.00794, 'HB1' : 1.00794, 'HB2' : 1.00794, 'HB3' : 1.00794,
               'HD1' : 1.00794, 'HD11': 1.00794, 'HD12': 1.00794, 'HD13': 1.00794, 'HD2' : 1.00794, 'HD21': 1.00794, 'HD22': 1.00794, 'HD23': 1.00794, 'HD3' : 1.00794,
               'HE'  : 1.00794, 'HE1' : 1.00794, 'HE2' : 1.00794, 'HE21': 1.00794, 'HE22': 1.00794, 'HE3' : 1.00794,
               'HG'  : 1.00794, 'HG1' : 1.00794, 'HG11': 1.00794, 'HG12': 1.00794, 'HG13': 1.00794, 'HG2' : 1.00794, 'HG21': 1.00794, 'HG22': 1.00794, 'HG23': 1.00794, 'HG3' : 1.00794,
               'HH'  : 1.00794, 'HH11': 1.00794, 'HH12': 1.00794, 'HH2' : 1.00794, 'HH21': 1.00794, 'HH22': 1.00794, 'HZ'  : 1.00794, 'HZ1' : 1.00794, 'HZ2' : 1.00794, 'HZ3' : 1.00794
             }


def help():
   line = textwrap.dedent("""\
                                   ##################################################
                                   ## Distance Fluctuation Analysis v 2.1 (Python) ##
                                   ##################################################
                                   Program developed by Giorgio Colombo Group
                             ---------------------------------------------------------------------------------------

                                         Copyright (c) 2022, Giorgio Colombo
		
                                   This work is licensed under a Apache 2.0 Licence

                             ---------------------------------------------------------------------------------------

                                   This script calculates the average distance fluctuation between residues in a MD simulation.
                                   It then returns the aminoacids showing a long-range coordinated motion.
                                   It is possible to use as input a .xyz (only C-alpha) file, or a .pdb file (only C-alpha,
                                   no-hydrogen or full atoms).
                                   Only the protein must be present in the structure files (no ligands, cofactors, waters or ions
                                   are allowed).
                                   It allows to calculate the distance fluctuations for the C-alphas or the sidechains.
                                   In the case of the sidechains distance fluctuation analysis it uses the center of mass of the
                                   residue sidechain, comprised of:
                                      C-alpha, other heavy atoms and eventually hydrogens (if present).
                               
                                   N.B. If the analysis is conducted only on C-alphas a better performance is obtained using either
                                        a. xyz file or a .pdb file with only c-alphas (and the specific flag '-p c-a'.
                                   
                                   N.B. The run time of the program is faster for .pdb files without hydrogens.

                                   N.B. The program requires python version 3.0 or higher to run

                                   Usage:
                                            python distance_fluctuation.py
                               """)
   return line

def valid_file(param):
    base, ext = os.path.splitext(param)
    if ext.lower() not in ('.pdb'):
        raise argparse.ArgumentTypeError('\n\nError, file must be a .pdb trajectory\n\n')
        sys.exit()
    return param

def get_args():
   parser=argparse.ArgumentParser(prog="distance_fluctuation.py",
                                  formatter_class=argparse.RawDescriptionHelpFormatter,
                                  description=help())
   required = parser.add_argument_group('Required arguments')
   optional = parser.add_argument_group('Optional arguments')
   required.add_argument('-ext', '--file_ext',  help = 'Input trajectory file iformat (options: .pdb or .xyz)',         choices = ['pdb', 'xyz'],                                                   type=valid_file)
   required.add_argument('-n',   '--n_aa',      help = 'Number of aminoacids in each frame of the .pdb trajectory',                                                                                 type = int     )
   optional.add_argument('-i',   '--in_name',   help = 'Name of the trajectory file (default: trj)',                    default = "trj",                                                            type = str     )
   optional.add_argument('-s',   '--s_frame',   help = 'Index of the initial frame (default = 0)',                      default = 0,                                                                type = int     )
   optional.add_argument('-e',   '--e_frame',   help = 'Index of the last frame (default = end)',                       default = 0,                                                                type = int     )
   optional.add_argument('-c',   '--cutoff',    help = 'Cutoff for the distance fluctuation analysis (default = 5)',    default = 5,                                                                type = int     )
   optional.add_argument('-t',   '--tolerance', help = 'Tolerance for the distance fluctuation analysis (default = 0)', default = 0,                                                                type = float   )
   optional.add_argument('-p',   '--pdb_type',  help = 'Type of pdb file submitted (only C-alpha, all protein atoms)',  default = 'c-a',                     choices = ['all', 'c-a'],              type = str     )
   optional.add_argument('-l',   '--local',     help = 'Type of local cutoff used, sequence or volume (default = seq)', default = 's',                       choices = ['s', 'seq', 'v', 'volume'], type = str     )
   optional.add_argument('-r',   '--radius',    help = 'Radius of the area cutoff (default = 7A)',                      default = 7,                                                                type = float   )
   optional.add_argument('-rs',  '--res_start', help = 'Starting residue if the PDB does not start with res. num. 1',   default = 1,                                                                type = int     )
   optional.add_argument('-ro',  '--rms_out',   help = 'RMS distance output filename (without extension)',              default = 'rmsdist_out',                                                    type = str     )
   optional.add_argument('-ao',  '--avg_out',   help = 'Average distance output filename (without extension)',          default = 'avgdist_out',                                                    type = str     )
   optional.add_argument('-so',  '--seq_out',   help = 'Profile sequence output filename (without extension)',          default = 'profile_sequence',                                               type = str     )
   optional.add_argument('-da',  '--dist_an',   help = 'Distance analysis output filename (without extension)',         default = 'distance_analysis_out',                                          type = str     )
   optional.add_argument('--blocks',            help = 'File with domains borders', default="", type=str)
   optional.add_argument('-rb',  '--rms_b_out', help = 'RMS distance blocks output filename', default='rmsdist_out', type=str )
   args = parser.parse_args()
   parser.print_help()
   infile = args.file_ext
   in_name = args.in_name
   num = args.n_aa
   start_frame = args.s_frame
   end_frame = args.e_frame
   cutoff = args.cutoff
   tolerance = args.tolerance
   rms_out = args.rms_out
   avg_out = args.avg_out
   seq_out = args.seq_out
   dist_an = args.dist_an
   local = args.local
   radius = args.radius
   res_num = args.res_start
   block_filename = args.blocks
   rms_b_out = args.rms_b_out
   seq = False
   pdb_all = False
   if re.search('pdb', args.file_ext) and args.pdb_type is None:
      parser.error("--file_ext pdb requires the pdb type.\n  Options are:\n      - C-a for a C-alpha only .pdb\n      - all for all protein atom .pdb file")
   if re.search('pdb', args.file_ext) and re.search('all', args.pdb_type):
      pdb_all = True
   if re.search('s|seq', local):
      seq = True
   return infile, in_name, num, start_frame, end_frame, cutoff, tolerance, rms_out, avg_out, seq_out, dist_an, seq, pdb_all, radius, res_num, block_filename, rms_b_out


def cut_table(table, start_frame, end_frame, real_numb,num):
   if start_frame == 0 and end_frame == 0:
      end_frame = real_numb/num
      if real_numb%(num) != 0:
         print("Premature end of file, check your input and try again!")
         sys.exit()
   elif start_frame == 0 and end_frame != 0:
      table = table[: end_frame]
   elif start_frame != 0 and end_frame == 0:
      end_frame = real_numb/num
      if real_numb%(num) != 0:
         print("Premature end of file, check your input and try again!")
         sys.exit()
      table = table[start_frame:]
   elif start_frame != 0 and end_frame != 0:
      table = table[start_frame:end_frame]
   return table, end_frame


def read_pdb_CA_file(infile, num, start_frame, end_frame):
   frames = []
   print("\nProcessing the .pdb file" )
   table = []
   with open(infile, 'r') as pdb_file:
      my_line = np.zeros((num, 3))
      lines = pdb_file.readlines()
      filt_coords = []
      count = 0
      real_numb = 0
      for line in lines:
         if 'ATOM' in line:
             #STEFANO MODIFIES, AS DOESN'T WORK IN PDBs WHEREIN XYZ ARE < -99.99999
             #coords = line[30:56].split()
             #filt_coords.append([float(coords[0]), float(coords[1]), float(coords[2])])
             filt_coords.append([float(line[30:38]),float(line[38:46]),float(line[46:54])])
             count += 1
             real_numb += 1
             if num == count:
                table.append(filt_coords)
                count = 0
                filt_coords = []
   table, end_frame = cut_table(table, start_frame, end_frame, real_numb,num)
   return table, end_frame


def read_pdb_all_file(infile, num, start_frame, end_frame, res_num):
    # Function to read the file and divide it into frames.
    # It returns a vector containing the coordinates in each frame and a vector containing the names of the residues.
   print("\nProcessing the .pdb file" )
   table = []
   tot_atoms = 0
   line_numb = 0
   real_numb = 0
   frame_numb = 0
   if num > 1:
      num = num + res_num -1
   with open(infile) as pdb_file:
      resid = res_num
      current_res = 1
      my_line = np.zeros((num, 3))
      residue_coord = []
      res_atom_mass = []
      frame = []
      loop = 0
      for idx, line in enumerate(pdb_file):
         if 'CRYST1' not in line and 'END' not in line and 'MODEL' not in line and 'TER' not in line and 'ENDMDL' not in line:
            first_half_line = line[:22]
            second_half_line = line[22:]
            line_split1 = first_half_line.split()
            line_split2 = second_half_line.split()
            atom_type = line_split1[2]
            current_res = int(line_split2[0])
            coords = [float(line_split2[1]), float(line_split2[2]), float(line_split2[3])]
            if resid == current_res:
               res_atom_mass.append(masses_pdb[atom_type])
               residue_coord.append(coords)
            else:
               center_of_mass = np.average(residue_coord, axis = 0, weights=res_atom_mass)
               frame.append(center_of_mass)
               resid += 1
               residue_coord = []
               res_atom_mass = []
               first_half_line = line[:22]
               second_half_line = line[22:]
               line_split1 = first_half_line.split()
               line_split2 = second_half_line.split()
               atom_type = line_split1[2]
               current_res = int(line_split2[0])
               coords = [float(line_split2[1]), float(line_split2[2]), float(line_split2[3])]
               res_atom_mass.append(masses_pdb[atom_type])
               residue_coord.append(coords)
               real_numb += 1
         elif 'TER' in line:
            if current_res < num:
               resid += 1
               center_of_mass = np.average(residue_coord, axis = 0, weights=res_atom_mass)
               frame.append(center_of_mass)
               residue_coord = []
               res_atom_mass = []
               continue
         else:
            if idx != 0:
               if end_frame != 0:
                  frame_numb += 1
                  if resid >= num:
                     center_of_mass = np.average(residue_coord, axis = 0, weights=res_atom_mass)
                     frame.append(center_of_mass)
                  table.append(frame)
                  residue_coord = []
                  res_atom_mass = []
                  frame = []
                  resid = res_num
                  current_res = 1
                  if frame_numb ==end_frame:
                     break
               else:
                  if resid >= num:
                     center_of_mass = np.average(residue_coord, axis = 0, weights=res_atom_mass)
                     frame.append(center_of_mass)
                  frame_numb += 1
                  if len(frame) > 0:
                     table.append(frame)
                  residue_coord = []
                  res_atom_mass = []
                  frame = []
                  resid = res_num
                  current_res = 1
   if start_frame != 0 and end_frame != 0:
        print()
   elif start_frame != 0 and end_frame == 0:
        end_frame = frame_numb
   elif end_frame == 0 and end_frame != 0:
        print()
   else:
        end_frame = frame_numb
   #table, end_frame = cut_table(table, start_frame, end_frame, real_numb)
   print("  Done!\n")
   return table, end_frame


def read_xyz_file(infile, num, start_frame, end_frame):
    # Function to read the file and divide it into frames.
    # It returns a vector containing the coordinates in each frame and a vector containing the names of the residues.
   frames = []
   print("\nProcessing the .xyz file" )
   table = []
   tot_atoms = 0
   #chunk_count = 0
   filtered_lines = []
   with open(infile) as xyz_file:
      lines = xyz_file.readlines()
      for line in lines:
         if line.strip().isdigit():
            continue
         else:
            filtered_lines.append(line)

   for i in range(0, len(filtered_lines), num):
      chunk = filtered_lines[i:i + num]
      my_line = np.zeros((num, 3))
      for atom, line in enumerate(chunk):
        tot_atoms += 1
        coords = line.split()
        for i in range(1,4):
           my_line[atom, i - 1] = coords[i]
      table.append(my_line)
   if start_frame != 0 and end_frame != 0:
     table = table[start_frame : end_frame]
   elif start_frame != 0:
     table = table[start_frame:]
   elif end_frame != 0:
     table = table[:end_frame]
   if tot_atoms % num != 0:
      print("premature end of file")
      sys.exit()
   if end_frame == 0:
      end_frame = len(table) - start_frame#chunk_count
   print("  Done!\n")
   return table, end_frame


def calc_average_dist(frames, num):
   print("\nCalculating average distance matrix:")
   ave_dist = np.zeros((num, num))
   nframes = 0
   temp_dist = np.zeros((num, num))
   for idx, frame in enumerate(frames):
      my_dist = distance.cdist(frame, frame, 'euclidean')
      temp_dist = np.add(temp_dist, my_dist)
      nframes +=1
   ave_dist = np.divide(temp_dist, nframes)
   print("   Done!\n")
   return ave_dist, nframes

def build_logical_matrix_df(num):
   logical_matrix_df = np.ones((num, num))   # Matrix to filter the residues in the range 1-x, 1+x
   for i in range(num):
      for j in range(num):
         if abs(i-j) < 2:
            logical_matrix_df[i, j]=0
   return logical_matrix_df

def calc_nearcutoff_seq(dist_fluct):
   print("\nCalculating the near cutoff:")
   nearcutoff = 0                             # Variable that define the nearcutoff
   counter = 0                                # Counter to calculate the nearcutoff
   # Compute the value of the nearcutoff
   for idx1, line in enumerate(dist_fluct):
      for idx2, elem in enumerate(line):
         if 0 < abs(idx1-idx2) < 5:
            nearcutoff += dist_fluct[idx1][idx2]
            counter += 1
   nearcutoff = nearcutoff / counter
   print("   Done.")
   print("      Nearcutoff value is: " + str(nearcutoff) + " A\n")
   return nearcutoff

def calc_nearcutoff_3d(average_dist, dist_fluct, radius, num):
   print("\nCalculating the near cutoff:")
   nearcutoff = 0
   counter = 0
   filt_avg = np.copy(average_dist)
   filt_avg[filt_avg > radius] = 0
   filt_avg[filt_avg > 0] = 1
   nc_matrix = np.multiply(dist_fluct, filt_avg)
   counter = np.sum(filt_avg)
   nearcutoff = np.sum(nc_matrix)
   nearcutoff_sum2 = np.sum(np.square(nc_matrix))
   nearcutoff_std = (nearcutoff_sum2/counter)-(nearcutoff/counter)**2
   nearcutoff = nearcutoff/counter
   print("   Done.")
   print("      Nearcutoff value is: " + str(nearcutoff) + " +/- " + str(nearcutoff_std) +  " A\n")
   return nearcutoff
   


def calc_local_fluct_seq(dist_fluct, num):
   # Calculate the matrix of the local fluctuations. N.B. To be changed for the 3D local fluctuations
   print("Calculating local sequence distance fluctuation: ", end="")
   logical_matrix_ps = np.zeros((num, num))   # Matrix to filter the residues further away than 3 from residue x  (Used to calculate the profile sequence distance fluctuation)
   local_fluctuation = np.zeros(num)          # Vector that will store the values for the sequence distance fluctuation
   # Initialize logical_matrix_nc
   for i in range(3, num-3):
      index = [i-2,i-1,i+1,i+2]
      for val in index:
         logical_matrix_ps[i, val]=1
   local_fluct = np.multiply(logical_matrix_ps, dist_fluct)
   local_fluctuation = np.sum(local_fluct, axis = 1)
   local_fluctuation = np.divide(local_fluctuation, 4.0)
   print("Done.")
   return local_fluctuation


def calc_dist_fluct(frames, num, ave_dist, nframes):
   print("Calculating distance fluctuation matrix:")
   dist_fluct = np.zeros((num, num))
   logical_matrix_df = build_logical_matrix_df(num)
   # Build and calculate the distance matrix 
   for frame in frames:
      my_dist = distance.cdist(frame, frame, 'euclidean')                           # Calculate the distance between the pairs of residues
      diff_matrix = np.subtract(my_dist, ave_dist)                                  # Subtract the average distance fluctuation to the local fluctuation)
      prod_matrix = np.square(diff_matrix)                                          # calculate the power of the difference between distance fluctuation and local fluctuation
      dist_fluct = np.add(dist_fluct, np.multiply(prod_matrix, logical_matrix_df))  # Compute the distance fluctuation matrix (diff_matrix * logical_matrix_df)

   dist_fluct = np.divide(dist_fluct, nframes)
   
   print("   Done!\n")
   return dist_fluct


def distance_analysis(ave_dist, dist_fluct, cutoff, nearcutoff, tolerance):
   print('Now doing the distance analysis')
   it = np.nditer(ave_dist, flags = ['multi_index'])
   ave_filt = ave_dist > cutoff
   ave_mask = ave_filt.astype(float)
   cut_tol = nearcutoff + tolerance
   
   dist_filt = dist_fluct < cut_tol
   dist_mask = dist_filt.astype(float)
   final_mask = np.multiply(ave_mask, dist_mask)
   coord_residues = np.multiply(dist_fluct, final_mask)
   print('   Done!\n')
   return coord_residues

def write_num(idx):
   num = ""
   if idx + 1 < 10:
      num = "   " + str(idx + 1)
   elif idx + 1 < 100:
      num = "  " + str(idx + 1)
   elif idx + 1 < 1000:
      num = " " + str(idx + 1)
   else:
      num = str(idx + 1)
   return num


def print_matrix_to_dat(in_matrix, out_filename, ndec, skip):
   size = in_matrix.shape
   with open(out_filename, 'w') as f:
      for i, line in enumerate(in_matrix):
         for j, elem in enumerate(line):
            if elem == 0 and skip:
               continue
            else:
               col = write_num(i)
               row = write_num(j)
               decimals = "." + str(ndec) + "f"
               val = format(elem, decimals)
               f.write(col + " " + row + "    " + val + "\n")
         f.write("\n")

def print_vector_to_dat(in_vect, out_filename, ndec, skip):
   size = len(in_vect)
   with open(out_filename, 'w') as f:
      for i in range(size):
         if in_vect[i] == 0 and skip:
            continue
         else:
            element = write_num(i)
            decimals = "." + str(ndec) + "f"
            val = format(in_vect[i], decimals)
            f.write(element + "     " + val + "\n")

def print_run_parameters(infile, in_name, num, start_frame, end_frame, cutoff, tolerance, rms_out, avg_out, seq_out, dist_an, seq, pdb_all, rms_b_out):
   print("\n\n")
   print("**--------------------------------------***--------------------------------------**")
   print("  Running the rmsdist script with the following parameters:")
   print("     - " + infile + ":                " + in_name)
   print("     - Residues (C-alpha): " + str(num))
   print("     - Start frame:        " + str(start_frame))
   if end_frame == 0:
      print("     - End frame:          Last")
   else:
      print("     - End frame:          " + str(end_frame))
   print("     - Cutoff:             " + str(cutoff))
   print("     - Tolerance:          " + str(tolerance))
   if seq:
      print("     - Local cutoff:       Sequence adjacency")
   else:
      print("     - Local cutoff:       3D adjacency")
   if pdb_all:
      print("     - Type of analysis:   Side chains")
   else:
      print("     - Type of analysis:   Main chain (C-alphas)")
#   print("     - Output files:       " + rms_out + ".dat\n                           " + avg_out +
#         ".dat\n                           " + seq_out + ".dat\n                           " + dist_an + ".dat")
   print("**--------------------------------------***--------------------------------------**")
   print()

def compute_blocks(filename, distance_fluctuation_matrix, num):
   print("Reducing the df matrix in blocks:")

   inf_borders = np.genfromtxt(filename,dtype=int,usecols=0)
   sup_borders = np.genfromtxt(filename,dtype=int,usecols=1)

   num_block = len(inf_borders)

   new_mat = np.zeros( (num_block,num_block) ,float)

   for i in range(num_block):
      for j in range(num_block):
         cnt = 0
         for x in range(inf_borders[i]-1,sup_borders[i],1):
            for y in range(inf_borders[j]-1,sup_borders[j],1):
               new_mat[i,j] += distance_fluctuation_matrix[x,y]
               cnt += 1
         new_mat[i,j] /= cnt

   print('   Done!\n')

   return new_mat


def main():

   infile, in_name, num, start_frame, end_frame, cutoff, tolerance, rms_out, avg_out, seq_out, dist_an, seq, pdb_all, radius, res_num, block_filename, rms_b_out= get_args() # Get the run variables for the program

   analysis_type = ""
   if re.search('pdb', infile):
      if not pdb_all:
         analysis_type = "_CA"
         in_name_noext = in_name.replace('.pdb', '')
         in_pdb_file = in_name_noext + ".pdb"
         print_run_parameters(infile, in_name, num, start_frame, end_frame, cutoff, tolerance, rms_out, avg_out, seq_out, dist_an, seq, pdb_all, rms_b_out)
         frames, end_frame = read_pdb_CA_file(in_pdb_file, num, start_frame, end_frame)
      else:
         analysis_type = "_sidechains"
         in_name_noext = in_name.replace('.pdb', '')
         in_pdb_file = in_name_noext + ".pdb"
         print_run_parameters(infile, in_name, num, start_frame, end_frame, cutoff, tolerance, rms_out, avg_out, seq_out, dist_an, seq, pdb_all, rms_b_out)
         frames, end_frame = read_pdb_all_file(in_pdb_file, num, start_frame, end_frame, res_num)

   elif re.search('xyz', infile):
      analysis_type = "_CA"
      in_name_noext = in_name.replace('.xyz', '')
      in_xyz_file = in_name_noext + ".xyz"
      print_run_parameters(infile, in_name, num, start_frame, end_frame, cutoff, tolerance, rms_out, avg_out, seq_out, dist_an, seq, pdb_all, rms_b_out)
      frames, end_frame = read_xyz_file(in_xyz_file, num, start_frame, end_frame)



   average_distance_matrix, nframes = calc_average_dist(frames, num)       # calculate the average distance matrix
   distance_fluctuation_matrix = calc_dist_fluct(frames, num, average_distance_matrix, nframes) # Calculate the distance fluctuation matrix and the nearcutoff
   local_fluctuation = calc_local_fluct_seq(distance_fluctuation_matrix, num)
   if seq:
      nearcutoff = calc_nearcutoff_seq(distance_fluctuation_matrix) # Nearcutoff calculated based on the sequence
   else:
      nearcutoff = calc_nearcutoff_3d(average_distance_matrix, distance_fluctuation_matrix, radius, num) # Change it to nearcutoff based on 3d sphere
   distance_analysis_matrix = distance_analysis(average_distance_matrix, distance_fluctuation_matrix, cutoff, nearcutoff, tolerance) # calculate the

   #Setting the names for the output dat files
   rms_out_filename = rms_out + "_" + str(start_frame) + "_" + format(end_frame, ".0f") + analysis_type + ".dat"
   avg_out_filename = avg_out + "_" + str(start_frame) + "_" + format(end_frame, ".0f") + analysis_type + ".dat"
   seq_out_filename = seq_out + "_" + str(start_frame) + "_" + format(end_frame, ".0f") + analysis_type + ".dat"
   dist_an_filename = dist_an + "_" + str(cutoff) + "_" + str(start_frame) + "_" + format(end_frame, ".0f") + analysis_type + ".dat"
   rms_b_out_filename = rms_b_out + "_b_" + str(cutoff) + "_" + str(start_frame) + "_" + format(end_frame, ".0f") + analysis_type + "_in_blocks_" + ".dat"

   #Check if bloks reduction is required
   if ( block_filename != ""):
      df_matrix_in_blocks = compute_blocks(block_filename, distance_fluctuation_matrix, num)


   #Saving the .dat files
   print("Saving output files.\n")
   print("  " + rms_out_filename, end = '', flush = True)
   print_matrix_to_dat(distance_fluctuation_matrix, rms_out_filename, 6, False)
   print(" ->  printed.")
   print("  " + avg_out_filename, end = '', flush = True)
   print_matrix_to_dat(average_distance_matrix,     avg_out_filename, 6, False)
   print(" ->  printed.")
   print("  " + seq_out_filename, end = '', flush = True)
   print_vector_to_dat(local_fluctuation,           seq_out_filename, 3, True)
   print(" ->  printed.")
   print("  " + dist_an_filename, end = '', flush = True)
   print_matrix_to_dat(distance_analysis_matrix,    dist_an_filename, 3, True)
   #np.savetxt("test.txt", distance_analysis_matrix)
   print("  " + rms_b_out_filename, end='', flush=True)
   print_matrix_to_dat(df_matrix_in_blocks, rms_b_out_filename, 3, True)

   print(" -> printed.\n\nAll done, enjoy!\n\n")


if __name__ == "__main__":
    main()
