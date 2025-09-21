#!/usr/bin/python3

##### Code by J. Nam
##### Last updated: 2025-07-17
##### Description: Over-representation analyses using parameters from config file

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
    ora_cfg = config["ora"]
    
    dir_out = ora_cfg["dir_out"]
    dir_log = ora_cfg["dir_log"]
    thread = ora_cfg["thread"]
    job_limit = ora_cfg["job_limit"]
    partition = ora_cfg["partition"]
    job_name = ora_cfg["job_name"]
    allgenes_table = ora_cfg["allgenes_table"]
    wgcna_obj_cancer = ora_cfg["wgcna_obj_cancer"]
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

    wrap_cmd = (
        f"time Rscript {project_root}/src/step20_ora.R "
        f"{dir_out} {allgenes_table} {wgcna_obj_cancer}"
    )
    
    subprocess.run(shlex.split(wrap_cmd), check=True)

if __name__ == "__main__":
    main()
