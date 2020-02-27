import ase
import ase.io
from ase import Atoms
import numpy as np
import os
import time
import multiprocessing as mp
# Ensuring correct version of spglib is imported
try:
    import spglib as spg
except ImportError:
    from pyspglib import spglib as spg
# Custom-made functions for workflow
#import symm_ops as so
#import visualisation_tools as vt
import config_equivalence as ce
import misc_tools as mt


def create_and_check_rand_async(i, td_atoms, oh_atoms, ox_atoms, orig_cell, orig_coords, symm_ops, symm_op_count, ase_cell_orig):
    """Workflow made into a function for compatibility with 'pool.apply_async'.
    Actions of workflow:
    - Creates random shuffles amongst lists of td and oh atoms in (i.e. substitutions in the alloy)
    - Uses this and data from the ase Atoms object of the original config to create a new ase Atoms object (rand_cfg)
    - Calls the function isEquiv to apply all symmetry operations of the parent config to the random structure, returns True if a match is found

    Args:
        i (int): Iterator from pool.apply_async, number of times to create and test a random config
        td_atoms (list): List of atomic numbers for atoms on the td sites in the original config
        oh_atoms (list): List of atomic numbers for atoms on the oh sites in the original config
        ox_atoms (list): List of atomic numbers of the oxygen atoms in the original config
        orig_cell (np array): Lattice vectors from ase Atoms object of original config (object.get_cell)
        orig_coords (np array): Fractional atomic positions from ase Atoms object of original config (object.get_scaled_positions())
        symm_ops (dictionary): Dictionary of translational and rotational symmetry operations of the parent config from spglib
        symm_op_count (int): Total number of symmetry operations for parent config
        ase_cell_orig (ase Atoms object): Ase Atoms object of the original config

    Returns:
        int: 0 if no matches are found 1 if any matches are found
    """
    degeneracy_count = 0
    # Randomly shuffle lists of td and oh atoms
    td_shuf = np.random.choice(td_atoms, size=td_atoms.shape)
    oh_shuf = np.random.choice(oh_atoms, size=oh_atoms.shape)
    all_atoms_shuf = np.concatenate((td_shuf, oh_shuf, ox_atoms), axis=0)
    # Create rand_cfg ase Atoms object with orig cfg coords and shuffled td and oh atomic numbers
    rand_cfg = Atoms(cell=orig_cell, scaled_positions=orig_coords, pbc=True)
    rand_cfg.set_atomic_numbers(all_atoms_shuf)
    isEquiv = ce.check_for_equiv(symm_ops, symm_op_count, ase_cell_orig, rand_cfg)
    if isEquiv:
        degeneracy_count += 1
    return degeneracy_count
def collect_result(result):
    """Collecting results from async parallel processors

    Args:
        result: takes result from 'pool.apply_async'

    Returns:
        Appends results from each sync process and outputs list of results to a global (defined in main as empty list)
    """
    global mp_results
    mp_results.append(result)


if __name__=='__main__':

    ### INPUTS:
    threshold = 1e-3 # Tolerance used by spglib to identify spacegroup
    scaling = 100 # scaling*total_combinations for random sampling of each config when searching for degeneracy
    data_locs = 'data/set_B.dat' # File where each line is location of all original (unrelaxed) POSCARs to be analysed
    struc_type = 'B' # 'A' when td sites fill first or 'B' when oh sites fill first
    inpt_file = 'data/setB_all.info' # Data file outputted from first processing step of workflow
    output_file = 'data/setB_all+degen.info' # New file to combine info from file above and that from this step of the worflow
    ### END OF INPUTS

    # Check available number of processes for parallelisation and let user set num to use
    print("Number of processors: ", mp.cpu_count())
    num_proc = int(input("Choose number of processors to use (<= number above): "))
    pool = mp.Pool(num_proc)
    mp_results = [] # Set global list to collect results from async 

    # Start the timer!
    t1 = time.time()
     
    with open(data_locs) as f:
        all_set_locs = f.readlines()      
    
    all_degen_counts = []
    for loc in all_set_locs:
        cfg_inpt = loc.rstrip()

        print('Analysing: '+cfg_inpt) #Debug    
        try:  
            ### Step 0: Read in orig config with ase (here it is an unrelaxed POSCAR from CASM that has been re-formatted to be readable by ase)
            orig_cfg = os.path.join(cfg_inpt, 'POSCAR_orig')
            ase_cell_orig = ase.io.read(orig_cfg, format='vasp')

            ### Step 1: Create parent (here choice to 'de-colour' original config so all TM's are Co )
            parent_cfg = mt.de_colour(ase_cell_orig, 'Co')

            ### Step 2: Obtain symmetry operations of parent
            spglib_cell = mt.get_spglib_from_ase(parent_cfg)
            symm_ops = spg.get_symmetry(spglib_cell, threshold)
            tot_symm_ops = [(r, t) for r, t in zip(symm_ops['rotations'], symm_ops['translations'])]
            symm_op_count = len(tot_symm_ops)

            ### Step 3a: Generate random substitutions of original config (respecting if cfg is A- or B-type)
            ### Step 3b: For each random cfg (one-at-a-time), apply all symm ops of parent and check for equivalence with orig cfg

            # Assigning count of Co_td and Co_oh based on if structure is set A or set B
            Co_count = ase_cell_orig.get_chemical_symbols().count('Co')
            # First check that config is not an end-member of the alloy
            if (Co_count == 0 or Co_count == 24):
                all_degen_counts.append(1)
                print('Symmetry degeneracy of alloy end-member is just 1.')
                continue # Move on to next config in data list, don't waste time with the rest of the analysis!
            else:
                if (struc_type == 'A'):
                    if (Co_count <= 8):
                        Co_td = Co_count
                        Co_oh = 0
                    elif (Co_count > 8):
                        Co_td = 8
                        Co_oh = Co_count - 8
                    else:
                        print('Error in checking Co_count')
                elif (struc_type == 'B'):
                    if (Co_count <= 16):
                        Co_oh = Co_count
                        Co_td = 0
                    elif (Co_count > 16):
                        Co_oh = 16
                        Co_td = Co_count - 16
                    else:
                        print('Error in checking Co_count')
                else:
                    print('Error checking if struc is set A or set B')

                # Generating random substitutions of orig cfg and checking for symm degeneracy with orig cfg
                orig_coords = ase_cell_orig.get_scaled_positions()
                orig_cell = ase_cell_orig.get_cell()
                all_atoms = ase_cell_orig.get_atomic_numbers()
                td_atoms = all_atoms[0:8]
                oh_atoms = all_atoms[8:24]
                ox_atoms = all_atoms[24:56]
                # Set attempts to be scaling*combination space (latter based on Co_td and Co_oh counts)
                # Intention of scaling is to increase likelihood that each possible substitution is sampled
                combinations = mt.calc_combs(Co_td, Co_oh)
                attempts = int(combinations*scaling)
                #attempts = 10 # reduce just for test phase
                #orig_atom_list = ce.atom_nums_with_coords_pdSorted(ase_cell) # Use when testing method with coord sorting
                # Timing comparison of configs run in parallel
                #t0 = time.time()
                for i in range(attempts):
                    pool.apply_async(create_and_check_rand_async, args=(i, td_atoms, oh_atoms, ox_atoms, orig_cell, orig_coords, symm_ops, symm_op_count, ase_cell_orig), callback=collect_result)
                #print('We compared '+str(attempts)+' configs')
                #print('It took {0} secs to compare cfgs'.format((time.time()-t0)))
                
                # Collect together results from all processors
                degeneracy_frac = float(sum(mp_results))/float(scaling) # Divide by scaling of combination space
        except:
            print('Error in processing config from: '+str(cfg_inpt))
        else: # If no errors in code above, add data to final degeneracy fractions list
            all_degen_counts.append(degeneracy_frac)      
            print('Scaled degeneracy count: '+str(degeneracy_frac)+', with: '+str(combinations)+' possible combinations.')

    pool.close() # Close after all files have been analysed
    pool.join() # Prevents next lines of code from being executed before all processors have finised

    ### Step 4: Add all_degen_counts list as extra column in .info files for setA or setB
    with open(inpt_file, 'r') as f_in:
        with open(output_file, 'w') as f_out:
            lines = f_in.readlines()
            f_out.write(lines[0].rstrip()+', symm_degen_frac\n')
            for line, symm_degen in zip(range(1, len(lines)), all_degen_counts):
                f_out.write(lines[line].rstrip()+' '+str(symm_degen)+'\n')
    
    print('')
    print('It took {0} secs to process the dataset'.format((time.time()-t1)))
