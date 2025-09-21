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
    predict_mirna_cfg = config["predict_mirna"]

    hubgenes_table = predict_mirna_cfg["hubgenes_table"]
    dir_out = predict_mirna_cfg["dir_out"]
    dir_log = predict_mirna_cfg["dir_log"]
    thread = predict_mirna_cfg["thread"]
    job_limit = predict_mirna_cfg["job_limit"]
    partition = predict_mirna_cfg["partition"]
    job_name = predict_mirna_cfg["job_name"]
    mir_lnc = predict_mirna_cfg["mir_lnc"]
    mir_gene1 = predict_mirna_cfg["mir_gene1"]
    mir_gene2 = predict_mirna_cfg["mir_gene2"]
    mir_gene3 = predict_mirna_cfg["mir_gene3"]
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

    wrap_cmd = (
        f"time Rscript {project_root}/src/step15_predict_mirna.R "
        f"{hubgenes_table} {dir_out} {mir_lnc} {mir_gene1} {mir_gene2} "
        f"{mir_gene3}"
    )
    
    subprocess.run(shlex.split(wrap_cmd), check=True)

if __name__ == "__main__":
    main()
