#!/usr/bin/python3

##### Code by J. Nam
##### Last updated: 2025-07-17
##### Description: Run CellRanger count using parameters from config file

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
    seurat_qc_cfg = config["quality_control"]

    seurat_obj = seurat_qc_cfg["seurat_obj"]
    dir_out = seurat_qc_cfg["dir_out"]
    dir_log = seurat_qc_cfg["dir_log"]
    thread = seurat_qc_cfg["thread"]
    job_limit = seurat_qc_cfg["job_limit"]
    partition = seurat_qc_cfg["partition"]
    job_name = seurat_qc_cfg["job_name"]
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

    wrap_cmd = (
        f"time Rscript {project_root}/src/step04_quality_control.R {seurat_obj} "
        f"{dir_out}"
    )
    
    subprocess.run(shlex.split(wrap_cmd), check=True)

if __name__ == "__main__":
    main()
