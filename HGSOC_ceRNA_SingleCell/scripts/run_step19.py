#!/usr/bin/python3

##### Code by J. Nam
##### Last updated: 2025-07-17
##### Description: Predicting shared miRNA targets using parameters from config file

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
    visualise_validation_cfg = config["visualise_validation"]

    dir_out = visualise_validation_cfg["dir_out"]
    dir_log = visualise_validation_cfg["dir_log"]
    thread = visualise_validation_cfg["thread"]
    job_limit = visualise_validation_cfg["job_limit"]
    partition = visualise_validation_cfg["partition"]
    job_name = visualise_validation_cfg["job_name"]
    fdr_sig = visualise_validation_cfg["fdr_sig"]

    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

    wrap_cmd = (
        f"time Rscript {project_root}/src/step19_visualise_validation.R "
        f"{dir_out} {fdr_sig}"
    )
    
    subprocess.run(shlex.split(wrap_cmd), check=True)

if __name__ == "__main__":
    main()
