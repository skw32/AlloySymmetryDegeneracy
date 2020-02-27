# Methods for creating plots to visualise atomic arrangements 

import matplotlib.pyplot as plt
import ase
import ase.io
from ase.visualize.plot import plot_atoms
from ase import Atoms
import numpy as np
from mpl_toolkits import mplot3d
from IPython.core.pylabtools import figsize
figsize(5, 5)


def plot_cell_ase(cell, title):
    """Use ase in-built function to create 2D plot
    NOTE: Not very nice for visualising a 3D structure and appears to always set cell vectors to positive values

    Args:
        cell (ase Atoms object): Ase atoms object for structure to visualise
        title (str): Set a title for the plot

    Returns:
        Produces figure of plot, originally used in development notebook with '%matplotlib inline'
    """
    fig, ax = plt.subplots()
    plt.title(title)
    plot_atoms(cell, ax, radii=0.3, rotation=('0x,0y,0z'))
      
      
def plot_cell_custom_3D(ase_cell, title):
    """Produces 3D plot of structure for a system
    NOTE: Can only produce a plot for a system with a maximum of 3 different element types

    Args:
        cell (ase Atoms object): Ase atoms object for structure to visualise
        title (str): Set a title for the plot

    Returns:
        Produces figure of plot, originally used in development notebook with '%matplotlib inline'
    """
    figsize(5, 5)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    lattice = ase_cell.get_cell()
    positions = ase_cell.get_positions()
    numbers = ase_cell.get_atomic_numbers()
    
    # Set axes limits as +/- max values of lattice vectors
    x_max = max(abs(lattice[0]))
    y_max = max(abs(lattice[1]))
    z_max = max(abs(lattice[2]))

    atom_types = np.unique(numbers)
    first_atom_type = atom_types[0]
    second_atom_type = atom_types[1]
    first_colour = 'b'
    second_colour = 'r'
    third_colour = 'g'
    for pos, species in zip(positions, numbers):
        x_coord, y_coord, z_coord = pos[0:3]
        atom_type = species
        if (atom_type == first_atom_type):
            colour = first_colour
        elif (atom_type == second_atom_type):
            colour = second_colour
        else:
            colour = third_colour
        x = 0.5 * np.outer(np.cos(u), np.sin(v)) + x_coord
        y = 0.5 * np.outer(np.sin(u), np.sin(v)) + y_coord
        z = 0.5 * np.outer(np.ones(np.size(u)), np.cos(v)) + z_coord
        ax.plot_surface(x, y, z, color=colour)
    #ax.set_xlim3d(-x_max, x_max)
    #ax.set_ylim3d(-y_max, y_max)
    #ax.set_zlim3d(-z_max ,z_max)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.title(title)
    plt.tight_layout()
    plt.show(block=False)