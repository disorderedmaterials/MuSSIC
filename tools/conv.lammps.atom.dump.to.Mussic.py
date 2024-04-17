import re
import math
import numpy as np
import multiprocessing as mp
import time
import sys
import argparse

def count_lines_streaming(filename):
  line_count = 0
  with open(filename, 'r') as f:
    for _ in f:  # Iterate over lines without storing them
      line_count += 1
  return line_count

def coord_lmp_to_xyz(lineStr, xLen, yLen, zLen):
    line = lineStr.split()
    xs, ys, zs = np.array(line[2:5], dtype=float) # scaled positions, floats
    return f"{line[1]} {xs * xLen:.8f} {ys * yLen:.8f} {zs * zLen:.8f}\n"

def process_lammps_file(input_file, output_file, everyNstep, cond="default"):
    """
    Reads a lammps trajectory file, extracts lattice parameters,
    and writes the modified content to a new file for MuSSIC XYZ.

    Args:
        input_file (str): Path to the input xyz file.
        output_file (str): Path to the output file.
    """
    #with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    success = True

    input_nb_lines = count_lines_streaming(input_file)

    infile = open(input_file, 'r')
    outfile = open(output_file, 'w')

    lineNb = 0
    stepNb = 0
    skipLines = 0
    coordLines_to_read = 0
    lastId = 0
    num_atoms_this_step = 0
    step_begin = True
    
    lattice_pattern = 'ITEM: BOX BOUNDS pp pp pp'
    particles_pattern = 'ITEM: ATOMS id type xs ys zs'
    timestep_pattern = 'ITEM: TIMESTEP'
    num_atoms_pattern = 'ITEM: NUMBER OF ATOMS'
    lattice_lohi_pattern = re.compile(r'(\d+\.\d+e[+-]\d+)\s(\d+\.\d+e[+-]\d+)')

    while True:
        line = infile.readline()
        if not line: # end of file
            break
        lineNb += 1
        progress = int(lineNb / input_nb_lines * 100)
        print(f"Progress: {progress}%", end="\r")
        
        #print(f"the current lineNb is {lineNb}")
        
        if step_begin and timestep_pattern in line:
            #stepNb = int(infile.readline())
            line = infile.readline()
            lineNb += 1
            
            if stepNb % everyNstep != 0: # Its a step to skip
                skipLines = int(num_atoms_this_step) + 7 # add the other lines
                try:
                    [next(infile) for _ in range(skipLines)]
                except:
                    # reach end of file
                    break
                lineNb += skipLines
                stepNb += 1
                continue
            
            line = infile.readline()
            lineNb += 1
            if num_atoms_pattern not in line  : 
                print("error timestep not followed by the number of atoms")
                success = False
                break
            
            num_atoms_this_step = infile.readline()
            lineNb += 1

            line = infile.readline()
            lineNb += 1

            if lattice_pattern not in line: 
                print("error number of atoms not followed by the lattice")
                success = False
                break
            
            # read the three lines of box coordinates
            xlo, xhi = np.array(lattice_lohi_pattern.search(infile.readline()).groups(), dtype=float)
            ylo, yhi = np.array(lattice_lohi_pattern.search(infile.readline()).groups(), dtype=float)
            zlo, zhi = np.array(lattice_lohi_pattern.search(infile.readline()).groups(), dtype=float)

            # convert to lengths
            xLen = abs(xlo-xhi)
            yLen = abs(ylo-yhi)
            zLen = abs(zlo-zhi)

            lineNb += 3

            outfile.write(f'{num_atoms_this_step}') # Write the number of atoms at the beginning of each steps
            
            # different syntax possible
            if cond == "default":
                outfile.write(f'{stepNb}\n') # default xyz format, write step number
            elif cond == "NPT":
                outfile.write(f'{xLen:.8f} {yLen:.8f} {zLen:.8f}\n') # NPT write box dimensions
            elif cond == "NVT":
                outfile.write('\n') # NVT write blank line
            else:
                print('Wrong cond option {cond}, should be "default" or "NVT" or "NPT"')
                success = False
                break

            line = infile.readline() # should be ITEM: ATOMS id type xs ys zs
            lineNb += 1
            
            if particles_pattern not in line: 
                print("error not listing particles after box dimensions")
                success = False
                break

            step_begin = False
            coordLines_to_read = int(num_atoms_this_step) # inform next loops to read coordinates
            
            continue
        # End of header
        
        elif coordLines_to_read > 0:
            # read in block and convert all coordinates in parallel
            
            lineNb += coordLines_to_read - 1 # all those lines will be read
            args0 = (line, xLen, yLen, zLen)
            args = [(infile.readline(), xLen, yLen, zLen) for _ in range(coordLines_to_read-1)]
            args.insert(0, args0)

            with mp.Pool(numThreads) as pool:
                convLines = pool.starmap(coord_lmp_to_xyz, args)
                outfile.write("".join(convLines))
            
            step_begin = True
            stepNb += 1
            coordLines_to_read = -1
        
        else:
            print(f"error expected beginning of step or coords to read, instead got {line} at line Nb {lineNb-1}")
            success=False
            break

    
    # End of file conversion - While True
    infile.close()
    outfile.close()
    
    return success
# End of function process lammps

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Converts a lammps atom dump to xyz')
    parser.add_argument('-input', required=True, type=str, help='Input dump file to convert')
    parser.add_argument('-every', required=False, type=int, default=1, help='Convert one step every this many steps')
    parser.add_argument('-cond', required=False, type=str, default="default", 
        choices={"default", "NVT", "NPT"}, help='XYZ type: default, NVT (for MuSSIC) or NPT (for MuSSIC)')
    parser.add_argument('-out', required=False, type=str, default="notDef", help='Output filename. If not specified, changes the extension to .formatted.xyz')
    parser.add_argument('-np', required=False, type=int, default=None, help='Number of parallel processes. Default is the number of cores available')
    args = parser.parse_args()

    input_file = args.input
    if args.out == "notDef":
        # If not specified convert to a files with appendix ".formatted"
        input_file_split = input_file.split(".")
        #extension = input_file_split[-1]
        filename = ".".join(input_file_split[:-1])
        output_file = f'{filename}.formatted.xyz'
    else:
        output_file = args.out
    everyNstep =  args.every
    cond = args.cond
    numThreads = args.np

    tic = time.perf_counter()
    result = process_lammps_file(input_file, output_file, everyNstep, cond)
    toc = time.perf_counter()
    if result:
        print(f"Converted successfully lammps file: {input_file} to: {output_file} in {toc - tic:0.4f} seconds")

