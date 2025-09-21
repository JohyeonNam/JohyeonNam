#!/usr/bin/python3

##### Code by J. Nam
##### Last updated: 2025-07-17
##### Description: Cell type annotation using parameters from config file

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
    seurat_annotation_cfg = config["singler_annotation"]

    seurat_obj = seurat_annotation_cfg["seurat_obj"]
    dir_out = seurat_annotation_cfg["dir_out"]
    dir_log = seurat_annotation_cfg["dir_log"]
    thread = seurat_annotation_cfg["thread"]
    job_limit = seurat_annotation_cfg["job_limit"]
    partition = seurat_annotation_cfg["partition"]
    job_name = seurat_annotation_cfg["job_name"]
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

    wrap_cmd = (
        f"time Rscript {project_root}/src/step10_singler_annotation.R "
        f"{seurat_obj} {dir_out}"
    )
    
    subprocess.run(shlex.split(wrap_cmd), check=True)

if __name__ == "__main__":
    main()
