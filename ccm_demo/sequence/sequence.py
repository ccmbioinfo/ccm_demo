import os
import shutil
import subprocess
import warnings

import torch
from Bio import Seq, SeqIO

from ccm_demo.sequence.utils import *

class Sequence:
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence
        self.device="cuda" if torch.cuda.is_available() else "cpu"


    def embeddings(self, model="esmc_300m", normalize=False):
        if model == "esmc_300m" or model == "esmc_g00m":
            embeddings=esm3_embeddings(sequence=self.sequence, model=model,
                                       normalize=normalize, device=self.device)
        else:
            raise NotImplementedError("That model is not implemented.")
        return embeddings


    def mutate(self, position, to, new_name=None):
        """
        position is 0 based
        """
        new_sequence = list(self.sequence)
        if position > len(new_sequence)-1:
            raise ValueError("Position is greater than the sequence length {}".format(len(new_sequence)))
        new_sequence[position] = to

        self.sequence = "".join(new_sequence)
        if new_name is not None:
            self.name = new_name
        return self

    def msa(self, database, destination, output_name="msa.a3m", cleanup=True):
        self.write(os.path.join(destination, "query.fasta"))
        script=os.path.join(__file__.replace("sequence.py", ""), "../scripts/run_mmseqs.sh")
        #TODO better packaging
        command=" ".join(["bash", script, " -i", str(os.path.abspath(os.path.join(destination, "query.fasta"))),
                            "-d", database, "-o", str(os.path.join(destination, output_name)),
                            "-c"])
        run=subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if run.returncode != 0:
            raise ChildProcessError("There was an error with mmseqs run see error below '\n' {}".format(run.stderr))
        else:
            print("mmseqs run complete the results are in {}/{}".format(destination, output_name))

        if cleanup:
            items=["query.fasta", "query.index", "query.dbtype", "query.source", "query.lookup", "query",
                   "query.tab", "query_h", "query_h.dbtype", "query_h.index", "align", "align.dbtype","align.index"]
            for item in items:
                os.remove(os.path.join(destination, item))


    def write(self, fpath):
        seq=SeqIO.SeqRecord(Seq.Seq(self.sequence), id=self.name, description="")
        with open(fpath, "w") as handle:
            SeqIO.write(seq, handle, "fasta")

    def __str__(self):
        return "Sequence with name {} and {} aas".format(self.name, len(self.sequence))

