#!/usr/bin/python3

##### Code by J. Nam
##### Last updated: 2025-07-17
##### Description: Harmony integration using parameters from config file

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
    seurat_harmony_cfg = config["seurat_harmony"]

    seurat_obj = seurat_harmony_cfg["seurat_obj"]
    dir_out = seurat_harmony_cfg["dir_out"]
    dir_log = seurat_harmony_cfg["dir_log"]
    thread = seurat_harmony_cfg["thread"]
    job_limit = seurat_harmony_cfg["job_limit"]
    partition = seurat_harmony_cfg["partition"]
    job_name = seurat_harmony_cfg["job_name"]
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

    wrap_cmd = (
        f"time Rscript {project_root}/src/step07_seurat_harmony.R "
        f"{seurat_obj} {dir_out}"
    )
    
    subprocess.run(shlex.split(wrap_cmd), check=True)

if __name__ == "__main__":
    main()
