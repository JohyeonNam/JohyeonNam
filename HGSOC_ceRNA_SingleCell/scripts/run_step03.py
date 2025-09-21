#!/usr/bin/python3

##### Code by J. Nam
##### Last updated: 2025-07-17
##### Description: Make seurat object from cellranger files using parameters from config file

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
    seurat_make_cfg = config["seurat_make"]

    dir_in = seurat_make_cfg["dir_in"]
    dir_out = seurat_make_cfg["dir_out"]
    dir_log = seurat_make_cfg["dir_log"]
    thread = seurat_make_cfg["thread"]
    job_limit = seurat_make_cfg["job_limit"]
    partition = seurat_make_cfg["partition"]
    job_name = seurat_make_cfg["job_name"]
    sample_sheet = seurat_make_cfg["sample_sheet"]
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

    wrap_cmd = (
        f"time Rscript {project_root}/src/step03_seurat_make.R {dir_in} "
        f"{dir_out} {sample_sheet}"
    )
    
    subprocess.run(shlex.split(wrap_cmd), check=True)

if __name__ == "__main__":
    main()
