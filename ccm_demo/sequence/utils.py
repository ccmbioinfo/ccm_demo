import torch

from esm.models.esmc import ESMC
from esm.sdk.api import ESMProtein, LogitsConfig

def esm3_embeddings(sequence, model, normalize=False, device="cuda"):
    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"

    if model not in ["esmc_300m", "esmc_600m"]:
        raise ValueError("Invalid model name")
    protein = ESMProtein(sequence)
    client = ESMC.from_pretrained(model).to(device)  # or "cpu"
    protein_tensor = client.encode(protein)
    logits_output = client.logits(
        protein_tensor, LogitsConfig(sequence=True, return_embeddings=True)
    )
    embeddings = logits_output.embeddings
    if normalize:
        embeddings = logits_output.embeddings[0].mean(dim=0)
    return embeddings

#TODO
def esm2_embeddings(model, normalize=False):
    pass