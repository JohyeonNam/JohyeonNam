#!/usr/bin/python3

##### Code by J. Nam
##### Last updated: 2025-07-17
##### Description: Project information to the full dataset using parameters from config file

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
    seurat_projection_cfg = config["seurat_projection"]

    seurat_obj = seurat_projection_cfg["seurat_obj"]
    dir_out = seurat_projection_cfg["dir_out"]
    dir_log = seurat_projection_cfg["dir_log"]
    thread = seurat_projection_cfg["thread"]
    job_limit = seurat_projection_cfg["job_limit"]
    partition = seurat_projection_cfg["partition"]
    job_name = seurat_projection_cfg["job_name"]
    dim = seurat_projection_cfg["dim"]
    res = seurat_projection_cfg["res"]
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

    wrap_cmd = (
        f"time Rscript {project_root}/src/step09_seurat_projection.R "
        f"{seurat_obj} {dir_out} {dim} {res}"
    )
    
    subprocess.run(shlex.split(wrap_cmd), check=True)

if __name__ == "__main__":
    main()
