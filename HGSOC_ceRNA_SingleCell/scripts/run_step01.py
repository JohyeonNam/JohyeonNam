#!/usr/bin/python3

##### Code by J. Nam
##### Last updated: 2025-07-17
##### Description: Run CellRanger mkref using parameters from config file

import yaml
from src.step01_cellranger_mkref import run_cellranger_mkref

def load_config(path):
    with open(path, "r") as f:
        return yaml.safe_load(f)

def main():
    config = load_config("config/params.yaml")
    cellranger_mkref_cfg = config["cellranger_mkref"]

    fasta_in = cellranger_mkref_cfg["fasta_in"]
    gtf_in = cellranger_mkref_cfg["gtf_in"]
    dir_out = cellranger_mkref_cfg["dir_out"]
    dir_log = cellranger_mkref_cfg["dir_log"]
    thread = cellranger_mkref_cfg.get("thread", 4)
    job_name = cellranger_mkref_cfg["job_name"]
    job_limit = cellranger_mkref_cfg["job_limit"]
    partition = cellranger_mkref_cfg["partition"]

    run_cellranger_mkref(
        fasta_in=fasta_in,
        gtf_in=gtf_in,
        dir_out=dir_out,
        dir_log=dir_log,
        thread=thread,
        job_name=job_name,
        fob_limit=job_limit,
        partition=partition
    )

if __name__ == "__main__":
    main()