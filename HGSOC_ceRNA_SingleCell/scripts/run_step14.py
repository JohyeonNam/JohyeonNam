#!/usr/bin/python3

##### Code by J. Nam
##### Last updated: 2025-07-17
##### Description: Co-expression analysis using parameters from config file

import subprocess
import yaml
import os
import shlex
from src.utils.mods_jnam import job_count_sleep, current_time

def load_config(path):
    with open(path, "r") as f:
        return yaml.safe_load(f)

def main():
    config = load_config("config/params.yaml")
    marker_hdwgcna_cfg = config["hdwgcna_coexpression"]

    seurat_obj = marker_hdwgcna_cfg["seurat_obj"]
    dir_out = marker_hdwgcna_cfg["dir_out"]
    dir_log = marker_hdwgcna_cfg["dir_log"]
    thread = marker_hdwgcna_cfg["thread"]
    job_limit = marker_hdwgcna_cfg["job_limit"]
    partition = marker_hdwgcna_cfg["partition"]
    job_name = marker_hdwgcna_cfg["job_name"]
    pcg_markers = marker_hdwgcna_cfg["pcg_markers"]
    lnc_markers = marker_hdwgcna_cfg["lnc_markers"]
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

    wrap_cmd = (
        f"time Rscript {project_root}/src/step14_hdwgcna_coexpression.R "
        f"{seurat_obj} {dir_out} {pcg_markers} {lnc_markers}"
    )
    
    subprocess.run(shlex.split(wrap_cmd), check=True)

if __name__ == "__main__":
    main()
