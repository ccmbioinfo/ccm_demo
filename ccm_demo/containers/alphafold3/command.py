# many thanks to @r-shadoff for writing most of this code

import json
import string
from ccm_demo.containers.utils import *

def generate_alpha_id(n):
    alpha_IDs = list(string.ascii_uppercase)
    if n > 26:  # Check if AA-ZZ should be appended to alpha_IDs
        for first in string.ascii_uppercase:
            for second in string.ascii_uppercase:
                alpha_IDs.append(first + second)
                if len(alpha_IDs) >= n:  # Ensures unnecessary IDs are not generated
                    return alpha_IDs
    return alpha_IDs[:n]

#TODO modifications
def generate_json(name, sequences, types, stoichiometry, fpath, seed=1, version=1):
    # to make it compatible with below
    if isinstance(stoichiometry, int):
        stoichiometry = [stoichiometry]

    if isinstance(types, str):
        types = [types]

    if isinstance(sequences, str):
        sequences = [sequences]

    total_chains = sum(map(int, stoichiometry))
    alpha_IDs = generate_alpha_id(total_chains)

    dialect = "alphafold3"
    if len(sequences) != len(stoichiometry):
        raise ValueError("sequences and stoichiometry must have same length")
    for item_type in types:
        if item_type not in ["protein", "ligand", "dna", "rna"]:
            raise ValueError("type must be one of 'protein', 'ligand', 'dna', 'rna'")

    protein_sequences = []
    alpha_index = 0

    for seq, item_type, stoich in zip(sequences, types, stoichiometry):
        ids = alpha_IDs[alpha_index:alpha_index + stoich]
        id_value = ids[0] if len(ids) == 1 else ids  # Convert to string if one item
        protein_sequences.append({item_type: {"id": id_value, "sequence": seq}})
        alpha_index += stoich

    json_file = {"name": name, "modelSeeds": [412], "sequences": protein_sequences,
                 "dialect": dialect, "version": version}
    file_name = f"{name}.json"
    json_path = os.path.abspath(os.path.join(fpath, file_name))
    with open(json_path, "w") as f:
        json.dump(json_file, f)

# this is for demo purposes
command_dict = {"bind_mounts": [{"local": "/hpf/largeprojects/ccmbio/acelik_files/alphafold3",
                                  "target": "/af3"},
                                 {"local":None,
                                  "target":"/inputs"},
                                 {"local":None,
                                  "target":"/outputs"}],
                "singularity_args":["--nv", "--cleanenv"],
                 "run_command":"""
                 python /af3/run_alphafold.py \\
                    --db_dir /af3/af3_db \\
                    --input_dir /inputs \\
                    --json_path /inputs/$file \\
                    --model_dir /af3/af3_weights \\
                    --output_dir /outputs \\
                    --run_data_pipeline={} \\
                    --run_inference={}
                 """}
