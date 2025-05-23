from ccm_benchmate.container_runner.container_runner import ContainerRunner
from esm.models.esm3 import ESM3
from esm.sdk.api import ESMProtein, SamplingConfig

def get_esm3_embeddings(file, normalize=True, device="cuda"):
    ESM3InferenceClient = ESM3.from_pretrained("esm3-open").to(device)

    protein = ESMProtein.from_pdb(file)
    protein_tensor = ESM3InferenceClient.encode(protein)

    if normalize:
        output = ESM3InferenceClient.forward_and_sample(
            protein_tensor, SamplingConfig(return_mean_embedding=True)
        )
        embeddings = output.mean_embedding
    else:
        output = ESM3InferenceClient.forward_and_sample(
            protein_tensor, SamplingConfig(return_per_residue_embeddings=True)
        )
        embeddings = output.per_residue_embedding
    return embeddings


def create_AF3_json(sequence_dict, use_msa=False):
    pass


def cleanup():
    pass

def setup_openmm():
    pass

def bio_emu_runner():
    pass




class Simulation:
    def __init__(self, protein):
        self.protein = protein

    def cleanup(self):
        pass

    # TODO openmm with a warning
    def simulate(self, use_bioemu=True, relax=False):
        pass

    def save(self):
        pass


