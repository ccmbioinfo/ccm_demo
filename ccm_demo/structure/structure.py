import os
import subprocess
import string
import json

from biotite.structure import io
import torch
import requests

from esm.models.esm3 import ESM3
from esm.sdk.api import ESMProtein, SamplingConfig


class Structure:
    def __init__(self, pdb):
        self.pdb = pdb
        self.sasa=None,
        self.embeddings=None
        self.foldseek=None
        self.device = "cuda" if torch.cuda.is_available() else "cpu"
        self.ESM3InferenceClient = ESM3.from_pretrained("esm3-open").to(self.device)

    def download(self, id, source="PDB", destination=None, load_after_download=True):
        if source == "PDB":
            url="http://files.rcsb.org/download/{}.pdb".format(id)
        elif source == "AFDB":
            url="https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v4.pdb".format(id)
        else:
            raise NotImplementedError("We can only download structures from PDB or AFDB")

        download = requests.get(url, stream=True)
        download.raise_for_status()
        with open("{}/{}.pdb".format(destination, id), "wb") as f:
            f.write(download.content)

        self.pdb="{}/{}.pdb".format(destination, id)

        if load_after_download:
            atom_array = io.load_structure(self.pdb)
        self.structure=atom_array
        return self

    def embeddings(self, model="esm3", normalize=False):
        protein = ESMProtein.from_pdb("./1utn.pdb")
        protein_tensor = self.ESM3InferenceClient.encode(protein)

        if normalize:
            output = self.ESM3InferenceClient.forward_and_sample(
                protein_tensor, SamplingConfig(return_mean_embedding=True)
            )
        else:
            output = self.ESM3InferenceClient.forward_and_sample(
                protein_tensor, SamplingConfig(return_per_residue_embeddings=True)
            )
        self.embeddings=output


    def align(self, other, destination):
        if self.pdb is None or other.pdb is None:
            raise ValueError("Cannot align structures without a PDB")

        command=["mustang", "-i", os.path.abspath(self.pdb), os.path.abspath(other.pdb), "-o", os.path.abspath(destination), "-r", "ON"]
        process=subprocess.run(command)
        if process.returncode != 0:
            raise ValueError("There was an error aligning structures. See error below \n {}".format(process.stderr))

        aligned_pdb=os.path.abspath(destination+".pdb")
        rotation_file=os.path.abspath(destination+".rot")
        return aligned_pdb, rotation_file

    def predict():
    """
    """
    # Structure Details - Change these for non-demo purposes
    structure_details={
        "complex_name":"pyr1-dimer_arath",
        "sequence":"MPSELTPEERSELKNSIAEFHTYQLDPGSCSSLHAQRIHAPPELVWSIVRRFDKPQTYKHFIKSCSVEQNFEMRVGCTRDVIVISGLPANTSTERLDILDDERRVTGFSIIGGEHRLTNYKSVTTVHRFEKENRIWTVVLESYVVDMPEGNSEDDTRMFADTVVKLNLQKLATVAEAMARNSGDGSGSQVT",
        "stoichiometry":1
    }
    complex_name = structure_details["complex_name"]
    sequence = structure_details['sequence']
    stoichiometry = structure_details["stoichiometry"]
   
   # All Required Directories - Change these for non-demo purposes
    required_dirs={
        # AF3 Directories
        "af3_dir":"/hpf/largeprojects/ccmbio/acelik_files/alphafold3",
        "af3_sif":"/hpf/largeprojects/ccmbio/acelik_files/alphafold3/alphafold3.sif",
        "af3_script":"/hpf/largeprojects/ccmbio/acelik_files/alphafold3/run_alphafold.py",
        "db_dir":"/hpf/largeprojects/ccmbio/acelik_files/alphafold3/af3_db",
        "model_dir":"/hpf/largeprojects/ccmbio/acelik_files/alphafold3/af3_weights",
        # Input/Output Directories
        "input_dir": "/hpf/largeprojects/ccmbio/rshadoff/ComplexPortalB/test",
        "output_dir": "/hpf/largeprojects/ccmbio/rshadoff/ComplexPortalB/test"
    }
    input_dir = required_dirs['input_dir']
    output_dir = required_dirs['output_dir']

    # Check inputs for validity 
    try:
        stoichiometry = int(stoichiometry)
    except ValueError:
        raise ValueError("Stoichiometry must be an integer.")
    
    for name, path in required_dirs.items():
        if path.endswith((".sif", ".py")):
            if not os.path.isfile(path):
                raise FileNotFoundError(f"{name} is not a valid file path: {path}")
        else:
            if not os.path.isdir(path):
                raise NotADirectoryError(f"{name} is not a valid directory: {path}")

        
    #Generate Alphabetical IDs 
    def generate_alpha_IDs(n):
        """"
        Generates a list of letters from A-Z, AA-ZZ based on the number of interacting proteins.
        :param n: The length of the list of interacting proteins
        """
        alpha_IDs = list(string.ascii_uppercase)  # list from A-Z generated automatically
        
        if n > 26:  # Check if AA-ZZ should be appended to alpha_IDs
            for first in string.ascii_uppercase:
                for second in string.ascii_uppercase:
                    alpha_IDs.append(first + second)
                    if len(alpha_IDs) >= n:  # Ensures unnecessary IDs are not generated
                        return alpha_IDs
        return alpha_IDs[:n]  # Returns the appropriate range of IDs
    alpha_IDs = generate_alpha_IDs(stoichiometry)
    ids = alpha_IDs[:stoichiometry]
    id_value = ids[0] if len(ids) == 1 else ids  # Convert to string if one item

    # Create JSON file for AF3
    json_file = {
        "name": complex_name,
        "modelSeeds": [412],
        "sequences": [
            {"protein": {
                "id": id_value,
                "sequence": sequence
            }
        }
        ],
        "dialect": "alphafold3",
        "version": 1
    }
    file_name = f"{complex_name}.json"
    json_path = os.path.join(input_dir, file_name) 
    with open(json_path, "w") as f:
        json.dump(json_file, f)
    
    # Create job script to run AF3
    job_script_path = os.path.join(input_dir, f"{complex_name}job.sh")
    with open(job_script_path, "w") as f:
        f.write(f"""#!/bin/bash

#SBATCH -N 1
#SBATCH -c 9
#SBATCH --mem 120G
#SBATCH --time=96:00:00
#SBATCH --tmp=800G
#SBATCH --output={job_script_path}.out
#SBATCH --error={job_script_path}.err

module load Singularity
files=$(ls {input_dir} | grep "json")

cp -R {required_dirs["af3_dir"]} $SCRATCH

for file in $files; do
if [[ -f {input_dir}/$file ]]; then
    echo "Processing JSON file: $file"

    singularity exec --cleanenv \\
        -B $SCRATCH/alphafold3/:/af3 \\
        -B {input_dir}:/inputs \\
        -B {output_dir}:/outputs \\
        {required_dirs["af3_sif"]} \\
        python /af3/run_alphafold.py \\
        --db_dir /af3/af3_db \\
        --input_dir /inputs \\
        --json_path /inputs/$file \\
        --model_dir /af3/af3_weights \\
        --output_dir /outputs \\
        --run_data_pipeline=False \\
        --run_inference=True
fi
done
""")

    try: 
        subprocess.run(["sbatch", job_script_path])
    except OSError:
        raise OSError(f"Unable to submit job with {job_script_path}")

class ProteinComplex:
    def __init__(self, sequences, stoichiometry, pdb=None):
        self.sequences = sequences
        self.stoichiometry = stoichiometry
        self.pdb = pdb

    def read_pdb(self, path):
        pass

    def contacts(self):
        pass

    def predict(self, use_msa=False, model="alphafold3"):
        # Structure Details - Change these for non-demo purposes
        structure_details={
            "complex_name": "ctx_crodu",
            "sequences": [
                "SSYGCYCGAGGQGWPQDASDRCCFEHDCCYAKLTGCDPT",
                "QEDGEIVCGEDDPCGTQICECDKAAAICFRNSMDT",
                "QFSPENCQGESQPC",
                "HLLQFNKMIKFETRKNAVPFYAFYGCYCGWGGQGRPKDATDRCCFVHDCCYGKLAKCNTKWDIYRYSLKSGYITCGKGTWCEEQICECDRVAAECLRRSLSTYKNGYMFYPDSRCRGPSETC"
        ],
            "stoichiometry": [1, 1, 1, 1]
        }
        complex_name = structure_details["complex_name"]
        sequences = structure_details['sequences']
        stoichiometry = structure_details["stoichiometry"]
    
    # All Required Directories - Change these for non-demo purposes
        required_dirs={
            # AF3 Directories
            "af3_dir":"/hpf/largeprojects/ccmbio/acelik_files/alphafold3",
            "af3_sif":"/hpf/largeprojects/ccmbio/acelik_files/alphafold3/alphafold3.sif",
            "af3_script":"/hpf/largeprojects/ccmbio/acelik_files/alphafold3/run_alphafold.py",
            "db_dir":"/hpf/largeprojects/ccmbio/acelik_files/alphafold3/af3_db",
            "model_dir":"/hpf/largeprojects/ccmbio/acelik_files/alphafold3/af3_weights",
            # Input/Output Directories
            "input_dir": "/hpf/largeprojects/ccmbio/rshadoff/ComplexPortalB/test",
            "output_dir": "/hpf/largeprojects/ccmbio/rshadoff/ComplexPortalB/test"
        }
        input_dir = required_dirs['input_dir']
        output_dir = required_dirs['output_dir']

        # Check inputs for validity 
        try:
            total_chains = sum(map(int, stoichiometry))
        except ValueError:
            raise ValueError("Stoichiometry must be a list of integers.")
        
        if len(stoichiometry) != len(sequences):
            raise ValueError("Missing stoichiometry value(s): Ensure each chain has a stoichiometry value.")
        
        for name, path in required_dirs.items():
            if path.endswith((".sif", ".py")):
                if not os.path.isfile(path):
                    raise FileNotFoundError(f"{name} is not a valid file path: {path}")
            else:
                if not os.path.isdir(path):
                    raise NotADirectoryError(f"{name} is not a valid directory: {path}")

        # Generate Alphabetical IDs 
        def generate_alpha_IDs(n):
            """"
            Generates a list of letters from A-Z, AA-ZZ based on the number of interacting proteins.
            :param n: The length of the list of interacting proteins
            """
            alpha_IDs = list(string.ascii_uppercase)  # list from A-Z generated automatically
            
            if n > 26:  # Check if AA-ZZ should be appended to alpha_IDs
                for first in string.ascii_uppercase:
                    for second in string.ascii_uppercase:
                        alpha_IDs.append(first + second)
                        if len(alpha_IDs) >= n:  # Ensures unnecessary IDs are not generated
                            return alpha_IDs
            return alpha_IDs[:n]  # Returns the appropriate range of IDs
        alpha_IDs = generate_alpha_IDs(total_chains)
        
        protein_sequences = []
        alpha_index = 0 

        for seq, stoich in zip(structure_details['sequences'], structure_details['stoichiometry']):
            ids = alpha_IDs[alpha_index:alpha_index + stoich]
            id_value = ids[0] if len(ids) == 1 else ids  # Convert to string if one item
            protein_sequences.append({
                "protein": {
                    "id": id_value,
                    "sequence": seq
                }
            })
            alpha_index += stoich
        
        # Create JSON file for AF3
        json_file = {
            "name": complex_name,
            "modelSeeds": [412],
            "sequences": protein_sequences,
            "dialect": model,
            "version": 1
        }
        file_name = f"{complex_name}.json"
        json_path = os.path.join(input_dir, file_name) 
        with open(json_path, "w") as f:
            json.dump(json_file, f)
        
        
        if use_msa=False:
            # Create job script to run AF3
            job_script_path = os.path.join(input_dir, f"{complex_name}job.sh")
            with open(job_script_path, "w") as f:
                f.write(f"""#!/bin/bash

        #SBATCH -N 1
        #SBATCH -c 10
        #SBATCH --mem 120G
        #SBATCH --time=96:00:00
        #SBATCH --tmp=800G
        #SBATCH --output={job_script_path}.out
        #SBATCH --error={job_script_path}.err

        module load Singularity
        files=$(ls {input_dir} | grep "json")

        cp -R {required_dirs["af3_dir"]} $SCRATCH

        for file in $files; do
        if [[ -f {input_dir}/$file ]]; then
            echo "Processing JSON file: $file"

            singularity exec --cleanenv \\
                -B $SCRATCH/alphafold3/:/af3 \\
                -B {input_dir}:/inputs \\
                -B {output_dir}:/outputs \\
                {required_dirs["af3_sif"]} \\
                python /af3/run_alphafold.py \\
                --db_dir /af3/af3_db \\
                --input_dir /inputs \\
                --json_path /inputs/$file \\
                --model_dir /af3/af3_weights \\
                --output_dir /outputs \\
                --run_data_pipeline=False \\
                --run_inference=True
        fi
        done
        """)
        else:
            # Create job script to run AF3
            job_script_path = os.path.join(input_dir, f"{complex_name}job.sh")
            with open(job_script_path, "w") as f:
                f.write(f"""#!/bin/bash

        #SBATCH -N 1
        #SBATCH -c 10
        #SBATCH --mem 120G
        #SBATCH --time=96:00:00
        #SBATCH --tmp=800G
        #SBATCH --output={job_script_path}.out
        #SBATCH --error={job_script_path}.err

        module load Singularity
        files=$(ls {input_dir} | grep "json")

        cp -R {required_dirs["af3_dir"]} $SCRATCH

        for file in $files; do
        if [[ -f {input_dir}/$file ]]; then
            echo "Processing JSON file: $file"

            singularity exec --cleanenv \\
                -B $SCRATCH/alphafold3/:/af3 \\
                -B {input_dir}:/inputs \\
                -B {output_dir}:/outputs \\
                {required_dirs["af3_sif"]} \\
                python /af3/run_alphafold.py \\
                --db_dir /af3/af3_db \\
                --input_dir /inputs \\
                --json_path /inputs/$file \\
                --model_dir /af3/af3_weights \\
                --output_dir /outputs \\
                --run_data_pipeline=True \\
                --run_inference=True
        fi
        done
        """)

        try: 
            subprocess.run(["sbatch", job_script_path])
        except OSError:
            raise OSError(f"Unable to submit job with {job_script_path}")