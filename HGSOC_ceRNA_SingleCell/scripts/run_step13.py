#!/usr/bin/python3

##### Code by J. Nam
##### Last updated: 2025-07-17
##### Description: Identification of the marker genes using parameters from config file

import subprocess
import yaml
import os
import shlex
from src.utils.mods_jnam import job_count_sleep, current_time, gtf_parser

def load_config(path):
    with open(path, "r") as f:
        return yaml.safe_load(f)

def main():
    config = load_config("config/params.yaml")
    marker_identification_cfg = config["marker_identification"]

    seurat_obj = marker_identification_cfg["seurat_obj"]
    dir_out = marker_identification_cfg["dir_out"]
    dir_log = marker_identification_cfg["dir_log"]
    thread = marker_identification_cfg["thread"]
    job_limit = marker_identification_cfg["job_limit"]
    partition = marker_identification_cfg["partition"]
    job_name = marker_identification_cfg["job_name"]
    gtf_type = marker_identification_cfg["gtf_type"]
    gtf_file = marker_identification_cfg["gtf_file"]
    gtf_parser(gtf_type, gtf_file, dir_out)
    idnametype = os.path.join(dir_out, f"custom_geneinfo_{gtf_type}.csv")
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..")) 

    wrap_cmd = (
        f"time Rscript {project_root}/src/step13_marker_identification.R "
        f"{seurat_obj} {dir_out} {idnametype}"
    )
    
    subprocess.run(shlex.split(wrap_cmd), check=True)

if __name__ == "__main__":
    main()
