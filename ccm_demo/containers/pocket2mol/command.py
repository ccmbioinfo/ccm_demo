import os
import yaml

from ccm_demo.containers.utils import *

config={
    "model":{"checkpoint":"./ckpt/pretrained_Pocket2Mol.pt"},
    "sample":{
        "seed":42,
        "mask_init":True,
        "num_samples":None,
        "beam_size":500,
        "max_steps":50,
        "threshold":{
            "focal_threshold":0.5,
            "pos_threshold":0.25,
            "element_threshold":0.3,
            "hasatom_threshold":0.6,
            "bond_threshold":0.4
        }
    }
}

#TODO this need to be updated majorly
def generate_yaml(dict, numsamples=10, max_steps=50, fpath=None):
    dict["sample"]["num_samples"] = numsamples
    dict["sample"]["max_steps"] = max_steps
    if fpath is not None:
        return dict
    else:
        yaml.dump(dict, os.path.abspath(os.path.join(fpath, "config.yaml")))


command_dict={"bind_mounts":[{"local":None, #this is pocket2mol git repo
                             "target":"/work"},
                             {"local":None, #this is where the pdb is
                              "target":"/inputs"},
                             {"local":None,
                              "target":"/outputs"},],
              "singularity_args":["--nv", "--cleanenv"],
              "run_command":"""
              python work/sample_for_pdb.py \\
                --pdb_path /inputs/{} \\
                --center {center} \\
                --bbox_size {size} \\
                --config {config} \\
                --outdir {outdir}               
              """
}


