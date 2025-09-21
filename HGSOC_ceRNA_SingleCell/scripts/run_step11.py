#!/usr/bin/python3

##### Code by J. Nam
##### Last updated: 2025-07-17
##### Description: CNV inference using parameters from config file

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
    infercnv_cfg = config["infercnv"]

    seurat_obj = infercnv_cfg["seurat_obj"]
    dir_out = infercnv_cfg["dir_out"]
    dir_log = infercnv_cfg["dir_log"]
    thread = infercnv_cfg["thread"]
    job_limit = infercnv_cfg["job_limit"]
    partition = infercnv_cfg["partition"]
    job_name = infercnv_cfg["job_name"]
    gene_order = infercnv_cfg["gene_order"]
    celltype_anno = infercnv_cfg["celltype_anno"]
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

    wrap_cmd = (
        f"time Rscript {project_root}/src/step11_infercnv.R "
        f"{seurat_obj} {dir_out} {gene_order} {celltype_anno}"
    )
    
    subprocess.run(shlex.split(wrap_cmd), check=True)

if __name__ == "__main__":
    main()
