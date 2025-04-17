# many thanks to @kobehuynh for writing most of this code

import os
import subprocess
import pandas as pd
from io import StringIO

from Bio.PDB import *
parser = PDBParser(PERMISSIVE=1)

from ccm_demo.binders.utils import *
from ccm_demo.containers.utils import *
from ccm_demo.containers.pocket2mol.command import *

from usearch_molecules.dataset import FingerprintedDataset, shape_ecfp4, shape_fcfp4, shape_maccs


class MoleculeBinder:
    def __init__(self, pdb):
        """
        the structure needed for the molecule binder
        :param protein: a Protein object, we need a structure, if the protein does not have a pdb file we can
        either use the structure module to predict one or download from pdb doesn't matter
        """
        self.pdb=pdb

    def find_pockets(self, **kwargs):
        """
        :param kwargs: these are additional key value pairs to be fed into fpocket, read its documentation for details
        :return:
        """
        command_params=[]
        for key, value in kwargs.items():
            command_params.append(f"--{key} {value}")

        command_params=" ".join(command_params)
        command="fpocket -f {pdb} -x -d {command_params}".format(pdb=self.pdb, command_params=command_params)
        run=subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if run.returncode != 0:
            raise RuntimeError(run.stderr.decode())
        else:
            results=run.stdout.decode()
            pocket_properties=pd.read_csv(StringIO(results), sep=" ")
            pocket_list=[file for file in os.listdir(self.pdb.replace(".pdb", "_out")) if file.endswith(".pdb")]
            pocket_coords=[get_pocket_dimensions(item) for item in pocket_list]
            return pocket_list, pocket_properties, pocket_coords

    #TODO implement fixed from KOBE
    def get_coord_cuboid(pdb_file, amino_acids=None, use_alpha_carbon=False):

        """
        Compute the coordinates of the 3D bounding cuboid for a set of amino acids in a PDB file

        Parameters:
        pdb_file: str
            Path to pdb file to be analyzed.

        amino_acids: list[str], optional, default=None
            List of amino acid residue names (3-letter codes, e.g, "ALA", "GLY") to be analyzed.
            If None, all residues are analyzed.

        use_alpha_carbon: boolean, default=False
            Whether to use only alpha carbons (Cα) for computing the bounding cuboid.
            If False, all atoms in the residue will be used.

        Output:
        Coordinates of the cuboid:
        -Xmax, Ymax, Zmax: float
        -Xmin, Ymin, Zmin: float

        Example:
        get_coord_cuboid("example.pdb", amino_acids=["ALA", "GLY"], use_alpha_carbon = True)
        """

        structure = parser.get_structure("protein", pdb_file)
        coord = []

        for model in structure:
            for chain in model:
                for residue in chain:
                    if amino_acids is None or residue.resname in amino_acids:
                        for atom in residue:
                            if use_alpha_carbon:
                                if atom.name == "CA":
                                    coord.append(atom.coord)
                            else:
                                coord.append(atom.coord)

        coord_numpy = np.array(coord)

        x_max, y_max, z_max = np.max(coord_numpy, axis=0)
        x_min, y_min, z_min = np.min(coord_numpy, axis=0)

        return {"xmax":x_max, "ymax":y_max, "zmax":z_max, "xmin":x_min, "ymin":y_min, "zmin":z_min}

    def generate_molecules(self, container, pocket, fpath, outdir, box_size=None, numsamples=10, max_steps=50):
        """
        this runs pocket2mol for a give pocket and generated n molecules. The pockets can be de
        :param pocket:
        :param n:
        :return:
        """
        generate_yaml(fpath, config=config, numsamples=numsamples, max_steps=max_steps)
        command_dict["bind_mounts"][0]["local"] = os.path.abspath("../containers/pocket2mol/Pocket2Mol")
        command_dict["bind_mounts"][1]["local"]=os.path.abspath(fpath)
        command_dict["bind_mounts"][2]["local"] = os.path.abspath(outdir)
        if box_size is None:
            box_size="23"
        command_dict["run_command"] = command_dict["run_command"].format(center=pocket, size=box_size,
                                                                         config="{}/config.yaml".format(fpath),
                                                                         outdir=outdir)
        runner=ContainerRunner(container, command_dict)
        runner.run()
        return runner


    def search_molecules(self, smiles, library, using="ecfp4", metric="tanimoto", n=100):
        """
        using useard molecules will search library for a given metric, the only metric implemented is tanimoto similarity
        :param library: a list of directories that contain usearch parquet files
        :param using: the implemented ones are ecfp4, fcfp4 and maccs
        :param metric: tanimoto otherwise NotImplementedError
        :param n: number of molecules to return
        :return: list of smiles for molecules
        """
        if metric !="tanimoto":
            raise NotImplementedError("metric must be tanimoto")

        if using not in ["ecfp4", "fcfp4", "maccs"]:
            raise NotImplementedError("method must be ecfp4 or fcfp4 or maccs")
        elif using=="ecfp4":
            shape=shape_ecfp4
        elif using=="fcfp4":
            shape=shape_ecfp4
        elif using=="maccs":
            shape=shape_maccs

        data=FingerprintedDataset(library, shape=shape)
        results=data.search(smiles=smiles, n)
        return results




    def calculate_properties(self, *args):
        """
        calculate numerous properties using RDkit simple wrapper
        :param args: properties to get
        :return: returns properites of a smiles string based on args
        """
        pass

    def verify(self, method="AF3"):
        """
        run AF3 to re-predict the structure to verify binding
        :param method:
        :param kwargs:
        :return:
        """
        pass


class PeptideBinder:
    def __init__(self, protein, region):
        self.protein = protein
        if region is not None:
            self.region = region

    def find_binding_site(self):
        """
        this will again use fpocket
        :return:
        """
        pass

    def generate_peptides(self, region, n):
        pass

    def verify(self, method="AF3", **kwargs):
        pass