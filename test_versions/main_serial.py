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
from symm_ops import *
from config_equivalence import *
from misc_tools import *
#from visualisation_tools import *


if __name__=='__main__':


    # Here just for a test_cfg in set A, later loop over whole of set in dir lists: 'set_A.dat' and 'set_B.dat'
    '''
    data_locs = 'set_A.dat'

    with open(data_locs) as f:
        all_set_locs = f.readlines()
        
    for loc in all_set_locs:
        inpt_cfg = loc.rstrip()
        print(inpt_cfg)
        
    '''

    # INPUTS:
    threshold = 1e-3 # used by spglib for identifying spacegroup of cfgs
    cfg_inpt = '/home/suzannekwallace/Projects/Co_xMn_{3-x}O_4/comparisonStudy_NNandMTP/oneElongAB+B1/block_B.2/con0.503'
    struc_type = 'B'
    inpt_file = 'data/setB_all.info'
    output_file = 'data/setB_all+degen.info'


    all_set_deg_frac = []
    try:    
        ### Step 0: Read in orig config with ase
        orig_cfg = os.path.join(cfg_inpt, 'POSCAR_orig')
        ase_cell_orig = ase.io.read(orig_cfg, format='vasp')

        ### Step 1: Create parent
        parent_cfg = de_colour(ase_cell_orig, 'Co')

        ### Step 2: Obtain symmetry operations of parent
        spglib_cell = get_spglib_from_ase(parent_cfg)
        symm_ops = spg.get_symmetry(spglib_cell, threshold)
        tot_symm_ops = [(r, t) for r, t in zip(symm_ops['rotations'], symm_ops['translations'])]
        symm_op_count = len(tot_symm_ops)

        ### Step 3a: Generate random substitutions of original config (respecting if cfg is A- or B-type)
        ### Step 3b: For each random cfg (one-at-a-time), apply all symm ops of parent and check for equivalence with orig cfg

        # Assigning count of Co_td and Co_oh based on if structure is set A or set B
        Co_count = ase_cell_orig.get_chemical_symbols().count('Co')
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
        degeneracy_count = 0
        orig_coords = ase_cell_orig.get_scaled_positions()
        orig_cell = ase_cell_orig.get_cell()
        all_atoms = ase_cell_orig.get_atomic_numbers()
        td_atoms = all_atoms[0:8]
        oh_atoms = all_atoms[8:24]
        ox_atoms = all_atoms[24:56]
        # Set attempts to be 100x combination space based on Co_td and Co_oh counts
        combinations = calc_combs(Co_td, Co_oh)
        attempts = int(combinations*100)
        #orig_atom_list = atom_nums_with_coords_pdSorted(ase_cell) # Use when testing method with coord sorting
        for i in range(attempts):
            counter = i
            # Randomly shuffle lists of td and oh atoms
            td_shuf = np.random.choice(td_atoms, size=td_atoms.shape)
            oh_shuf = np.random.choice(oh_atoms, size=oh_atoms.shape)
            all_atoms_shuf = np.concatenate((td_shuf, oh_shuf, ox_atoms), axis=0)
            # Create rand_cfg ase Atoms object with orig cfg coords and shuffled td and oh atomic numbers
            rand_cfg = Atoms(cell=orig_cell, scaled_positions=orig_coords, pbc=True)
            rand_cfg.set_atomic_numbers(all_atoms_shuf)
            # Checking timing for comparing cfgs
            #t0 = time.time()
            isEquiv = check_for_equiv(symm_ops, symm_op_count, ase_cell_orig, rand_cfg)
            #print('It took {0} secs to compare cfgs'.format((time.time()-t0)))
            if isEquiv:
                degeneracy_count += 1

        degeneracy_frac = float(degeneracy_count)/float(counter)


    except:
        print('Error in processing config from: '+str(cfg_inpt))
    else: # If no errors in code above, add data to final degeneracy fractions list
        all_set_deg_frac.append(degeneracy_frac)

        
    print('Degeneracy fraction: '+str(degeneracy_frac)+' of '+str(counter)+' attempts')


    '''
        ### Step 4: Add all_set_deg_frac list as extra column in .info files for setA or setB
        with open(inpt_file, 'r') as f_in:
            with open(output_file, 'w') as f_out:
                lines = f_in.readlines()
                f_out.write(lines[0].rstrip()+', symm_degen_frac\n')
                for line, symm_degen in zip(range(1, len(lines)), all_set_deg_frac):
                    f_out.write(lines[line].rstrip()+' '+str(symm_degen)+'\n')
                    
    '''     
