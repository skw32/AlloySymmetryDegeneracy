import ase
import ase.io
import math
from ase import Atoms
import numpy as np


# Methods for converting between ase and spglib cell formats
def get_spglib_from_ase(ase_cell):
    """Convert information from ase Atoms object to format read in by spglib

    Args:
        ase_cell (ase Atoms object): Structure to convert from ase Atoms object format

    Returns:
        np array: Input format of structure for spglib
    """
    lattice = ase_cell.get_cell()
    positions = ase_cell.get_scaled_positions()
    numbers = ase_cell.get_atomic_numbers()
    spglib_cell = (lattice, positions, numbers)
    return spglib_cell
def get_ase_from_spglib(spglib_cell):
    """Convert information from spglib format to ase Atoms object

    Args:
        spglib_cell (np array): Input format for spglib

    Returns:
        ase Atoms object: Input for spglib converted into ase Atom object
    """
    ase_cell = Atoms(cell=spglib_cell[0], scaled_positions=spglib_cell[1], pbc=True)
    ase_cell.set_atomic_numbers(spglib_cell[2])
    return ase_cell


# For preparing parent of given configuration
def de_colour(ase_cell, species):
    """Function used to convert all Mn into Co or all Co into Mn (i.e. 'de-colour' config to generate the parent for symmetry analysis)

    Args:
        ase_cell (ase Atoms object): Original config to make parent of
        species (str): Either 'Co' or 'Mn' depending on which species the alloy should be 'de-coloured' to

    Returns:
        ase Atoms object: Object modified so that all TMs are either Co or Mn
    """
    str_atoms = ase_cell.get_atomic_numbers().astype(str)
    if (species == 'Mn'):
        de_colour_atoms = np.char.replace(str_atoms, '27', '25')
    elif (species == 'Co'):
        de_colour_atoms = np.char.replace(str_atoms, '25', '27')
    else:
        print('Error in species selection, should be Co or Mn')
    int_atoms_de_colour = de_colour_atoms.astype(int)
    de_coloured_ase_cell = Atoms(cell=ase_cell.get_cell(), scaled_positions=ase_cell.get_scaled_positions(), pbc=True)
    de_coloured_ase_cell.set_atomic_numbers(int_atoms_de_colour)
    return de_coloured_ase_cell


def calc_combs(Co_td, Co_oh):
    """Calculate total number of combinations for fixed number of Co on td and oh sites
    n_td!/ r_td!*(n_td-r_td)! * n_oh!/ r_oh!*(n_oh-r_oh)!
    Where n is total number of td or oh sites and r is the number occupies by Co

    Args:
        Co_td (int): Number of Co on tetrahedral sites in the structure
        Co_oh (int): Number of Co on octahedral sites in the structure

    Returns:
        int: total number of combinations based on fixed number of Co on td and oh sites inputted
    """
    oh_combs = math.factorial(16)/(math.factorial(Co_oh)*math.factorial(16- Co_oh))
    td_combs = math.factorial(8)/(math.factorial(Co_td)*math.factorial(8-Co_td))
    total_combs = oh_combs*td_combs
    return total_combs