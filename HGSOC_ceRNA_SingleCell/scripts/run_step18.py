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
    tcga_coexpression_cfg = config["tcga_coexpression"]

    dir_out = tcga_coexpression_cfg["dir_out"]
    dir_log = tcga_coexpression_cfg["dir_log"]
    thread = tcga_coexpression_cfg["thread"]
    job_limit = tcga_coexpression_cfg["job_limit"]
    partition = tcga_coexpression_cfg["partition"]
    job_name = tcga_coexpression_cfg["job_name"]
    dir_tcga = tcga_coexpression_cfg["dir_tcga"]
    pcg_count = tcga_coexpression_cfg["pcg_count"]
    lnc_count = tcga_coexpression_cfg["lnc_count"]
    miRNA_count = tcga_coexpression_cfg["miRNA_count"]
    pred_cancer = tcga_coexpression_cfg["pred_cancer"]
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

    wrap_cmd = (
        f"time Rscript {project_root}/src/step18_tcga_coexpression.R "
        f"{dir_out} {dir_tcga} {pcg_count} {lnc_count} {miRNA_count} "
        f"{pred_cancer}"
    )
    
    subprocess.run(shlex.split(wrap_cmd), check=True)

if __name__ == "__main__":
    main()
