#!/usr/bin/python3

##### Code by J. Nam
##### Last updated: 2025-07-17
##### Description: Making TCGA expression count files using parameters from config file

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
    tcga_maketable_cfg = config["tcga_maketable"]

    dir_out = tcga_maketable_cfg["dir_out"]
    dir_log = tcga_maketable_cfg["dir_log"]
    thread = tcga_maketable_cfg["thread"]
    job_limit = tcga_maketable_cfg["job_limit"]
    partition = tcga_maketable_cfg["partition"]
    job_name = tcga_maketable_cfg["job_name"]
    dir_tcga = tcga_maketable_cfg["dir_tcga"]
    sample_sheet = tcga_maketable_cfg["sample_sheet"]
    clinical_sheet = tcga_maketable_cfg["clinical_sheet"]
    idnametype = tcga_maketable_cfg["idnametype"]
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

    wrap_cmd = (
        f"time Rscript {project_root}/src/step17_tcga_maketable.R "
        f"{dir_out} {dir_tcga} {sample_sheet} {clinical_sheet} {idnametype}"
    )
    
    subprocess.run(shlex.split(wrap_cmd), check=True)

if __name__ == "__main__":
    main()
