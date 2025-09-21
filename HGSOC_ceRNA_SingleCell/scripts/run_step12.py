#!/usr/bin/python3

##### Code by J. Nam
##### Last updated: 2025-07-17
##### Description: Pathway analyses using parameters from config file

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
    escape_cfg = config["escape"]

    seurat_obj = escape_cfg["seurat_obj"]
    dir_out = escape_cfg["dir_out"]
    dir_log = escape_cfg["dir_log"]
    thread = escape_cfg["thread"]
    job_limit = escape_cfg["job_limit"]
    partition = escape_cfg["partition"]
    job_name = escape_cfg["job_name"]
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

    wrap_cmd = (
        f"time Rscript {project_root}/src/step12_escape.R "
        f"{seurat_obj} {dir_out}"
    )
    
    subprocess.run(shlex.split(wrap_cmd), check=True)

if __name__ == "__main__":
    main()
