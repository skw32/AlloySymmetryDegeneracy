# Methods for applying symmetry operations to ase atoms objects

import ase
import ase.io
from ase import Atoms
import numpy as np


def rotate_cfg(ase_cell, rotate_symm_op):
    """Apply spglib rotation symmetry operation (matrix) to cell lattice vectors 

    Args:
        ase_cell (ase Atoms object): Cell to be transformed
        rotate_symm_op (np array): A single rotation matrix from the ['rotations'] dictionary key outputted by spglib

    Returns:
        ase Atoms object: The original ase Atoms object with the cell lattice vectors transformed
    """
    #rotated_cell = np.matmul(rotate_symm_op, ase_cell.get_cell())
    rotated_cell = np.matmul(ase_cell.get_cell(), rotate_symm_op) #test after reading spglib App. A, order doesn't matter here
    rotated = Atoms(cell=rotated_cell, scaled_positions=ase_cell.get_scaled_positions(), pbc=True)
    rotated.set_atomic_numbers(ase_cell.get_atomic_numbers())
    return rotated
def translate_cfg(ase_cell, translate_symm_op):
    """Apply the spglib translation symmetry operation (vector) to each atom

    Args:
        ase_cell (ase Atoms object): Cell to be transformed
        rotate_symm_op (np array): A single translation vector from the ['translations'] dictionary key outputted by spglib

    Returns:
        ase Atoms object: The original ase Atoms object with the atomic positions transformed
    """
    # Note: spglib uses fractional coords!
    positions = ase_cell.get_scaled_positions()
    translated_positions = positions + translate_symm_op
    translated = Atoms(cell=ase_cell.get_cell(), scaled_positions=translated_positions, pbc=True)
    translated.set_atomic_numbers(ase_cell.get_atomic_numbers())
    return translated

# A test function to account for periodicity in final translated and rotated structure (to ensure consistency with spglib)
def abs_latt_vecs(ase_cell):
    """Takes the absolute value of the lattice vectors to transform the structure back into the positive quadrant (account for periodicity)

    Args:
        ase_cell (ase Atoms object): Cell to be transformed

    Returns:
        ase Atoms object: Structure with an absolute value for the lattice vectors (i.e. no negative cell dimensions)
    """
    lattvec = ase_cell.get_cell()
    lattvec_abs = np.absolute(lattvec)
    ase_cell_abs = Atoms(cell=lattvec_abs, scaled_positions=ase_cell.get_scaled_positions(), pbc=True)
    ase_cell_abs.set_atomic_numbers(ase_cell.get_atomic_numbers())
    return ase_cell_abs

def all_operations(ase_cell, symm_ops, op_num):
    """Takes full symmetry operations dictionary outputted by spglib and applies both the translation and rotation of the same operation number

    Args:
        ase_cell (ase Atoms object): Cell to be transformed
        symm_ops (dictionary): Symmetry operations outputted by spglib with keys ['rotations'] and ['translations']
        op_num (int): Selected operation to perform from all possible operations in symm_ops dictionary

    Returns:
        ase Atoms object: Structure transformed by both symmetry operations of selected operation, taking into account periodicity
    """
    rotated_cell = rotate_cfg(ase_cell, symm_ops['rotations'][op_num])
    rotated_and_translated = translate_cfg(rotated_cell, symm_ops['translations'][op_num])
    rotated_and_translated_abs_cell = abs_latt_vecs(rotated_and_translated)
    return rotated_and_translated_abs_cell