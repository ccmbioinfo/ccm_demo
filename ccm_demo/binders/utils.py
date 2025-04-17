# many thanks to @kobehuynh for writing most of this code

import os
from Bio.PDB import PDBParser
import numpy as np

import ccm_demo.containers.alphafold3.command as af3 #this is needed for the verify commands


def get_pocket_dimensions(pocket_path):
    """
    Args:
        pocket_path (str): Path to the PDB file of the binding pocket.

    Returns:
        tuple:
            center (list of float): The [x, y, z] center coordinates.
            bbox_size (float): The maximum length in any dimension.
    """
    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure("pocket", pocket_path)

    coord = []

    # Extract all atom coordinates from the pocket structure
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    coord.append(atom.coord)

    # Convert list of coordinates to NumPy array
    coord_numpy = np.array(coord)

    # Find max and min along each axis
    x_max, y_max, z_max = np.max(coord_numpy, axis=0)
    x_min, y_min, z_min = np.min(coord_numpy, axis=0)

    # Compute the size of the bounding box (max extent)
    bbox_size = max(x_max - x_min, y_max - y_min, z_max - z_min)

    # Compute the geometric center of the pocket
    center = [
        (x_max + x_min) / 2,
        (y_max + y_min) / 2,
        (z_max + z_min) / 2
    ]
    return center, bbox_size

def generate_af3_runner():
    pass