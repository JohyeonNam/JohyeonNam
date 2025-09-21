#!/usr/bin/python3

##### Code by J. Nam
##### Last updated: 2025-07-17
##### Description: Run CellRanger count using parameters from config file

import yaml
from src.step02_cellranger_count import run_cellranger_count

def load_config(path):
    with open(path, "r") as f:
        return yaml.safe_load(f)

def main():
    config = load_config("config/params.yaml")
    cellranger_count_cfg = config["cellranger_count"]

    dir_in = cellranger_count_cfg["dir_in"]
    dir_out = cellranger_count_cfg["dir_out"]
    dir_log = cellranger_count_cfg["dir_log"]
    thread = cellranger_count_cfg["thread"]
    job_limit = cellranger_count_cfg["job_limit"]
    partition = cellranger_count_cfg["partition"]
    job_name = cellranger_count_cfg["job_name"]
    reference = cellranger_count_cfg["reference"]
    sample_sheet = cellranger_count_cfg["sample_sheet"]

    run_cellranger_count(
        dir_in=dir_in,
        dir_out=dir_out,
        dir_log=dir_log,
        thread=thread,
        job_limit=job_limit,
        partition=partition,
        job_name=job_name,
        reference=reference,
        sample_sheet=sample_sheet
    )

if __name__ == "__main__":
    main()
