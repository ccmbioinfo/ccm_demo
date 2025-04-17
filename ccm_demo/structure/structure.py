import os
import subprocess

import requests
import torch
from biotite.structure import io, to_sequence, superimpose, sasa
from biotite.structure.alphabet import to_3di

from ccm_demo.structure.utils import *
from ccm_demo.containers.alphafold3.command import *

class Structure:
    def __init__(self, pdb=None, sequence=None, predict=True, model="AF3"):
        self.pdb = pdb
        if pdb is not None:
            self.structure = io.load_structure(self.pdb)
        if sequence is not None and self.pdb is None:
            self.sequence = sequence
            if predict:
                self.structure = self.predict(self.sequence, model)
        if sequence is not None and pdb is not None:
            if self.sequence != to_sequence(self.structure):
                raise ValueError("Sequences do not match")
            else:
                self.sequence = sequence
        self.sasa = None,
        self.embeddings = None
        self.seq_3di = None
        self.device = "cuda" if torch.cuda.is_available() else "cpu"

    #TODO check for multiple chains
    def download(self, id, source="PDB", destination=None, load_after_download=True):
        if source == "PDB":
            url = "http://files.rcsb.org/download/{}.pdb".format(id)
        elif source == "AFDB":
            url = "https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v4.pdb".format(id)
        else:
            raise NotImplementedError("We can only download structures from PDB or AFDB")

        download = requests.get(url, stream=True)
        download.raise_for_status()
        with open("{}/{}.pdb".format(destination, id), "wb") as f:
            f.write(download.content)

        self.pdb = "{}/{}.pdb".format(destination, id)

        if load_after_download:
            atom_array = io.load_structure(self.pdb)
        self.structure = atom_array
        return self

    def calculate_embeddings(self, model="esm3", normalize=False):
        if model == "esm3":
            get_esm3_embeddings(self.pdb, normalize=normalize)
        else:
            raise NotImplementedError("We can only calculate embeddings from ESM3")

    def calculate_sasa(self, **kwargs):
        accessibility = sasa(self.structure, **kwargs)
        self.sasa = accessibility
        return self

    def get_3di(self):
        sequences, chain_starts = to_3di(self.structure)
        self.seq_3di = sequences
        return self

    #TODO replace with biotite
    def align(self, other, destination):
        if self.pdb is None or other.pdb is None:
            raise ValueError("Cannot align structures without a PDB")

        command = ["mustang", "-i", os.path.abspath(self.pdb), os.path.abspath(other.pdb), "-o",
                   os.path.abspath(destination), "-r", "ON"]
        process = subprocess.run(command)
        if process.returncode != 0:
            raise ValueError("There was an error aligning structures. See error below \n {}".format(process.stderr))

        aligned_pdb = os.path.abspath(destination + "results.pdb")
        rotation_file = os.path.abspath(destination + "results.rms_rot")
        html_report = os.path.abspath(destination + "results.html")
        return aligned_pdb, rotation_file, html_report

    #TODO implement AF2, and omegafold
    def predict(self, output_path, container, inference=True, pipeline=False, model="AF3"):
        if model != "AF3":
            raise NotImplementedError("We can only predict structures with AF3")
        generate_json(self.name, "protein", self.sequence, stoichiometry=1, fpath=output_path)
        os.makedirs(os.path.abspath(os.path.join(output_path, "af3_prediction")), exist_ok=True)
        command_dict["run_command"].format(pipeline, inference)
        command_dict["bind_mounts"][1] = os.path.abspath(output_path)
        command_dict["bind_mounts"][2]=os.path.abspath(os.path.join(output_path, "af3_prediction"))
        runner=ContainerRunner(container, command_dict)
        runner=runner.run()
        return runner

    def write(self, fpath):
        io.save_structure(self.pdb, fpath)


#TODO implement chai1 and boltz (later)
class ProteinComplex(Structure):
    def __init__(self, pdb=None, sequence=None):
        super().__init__(pdb=pdb, sequence=sequence)

    def get_chain(self, chain_id):
        pass

    def get_chain_ids(self):
        pass

    def contacts(self):
        pass


