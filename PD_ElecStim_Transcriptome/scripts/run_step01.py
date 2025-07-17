#!/usr/bin/python3

##### Code by J. Nam
##### Last updated: 2025-07-17
##### Description: Run TrimGalore using parameters from config file

import yaml
from src.step01_fastqc import run_fastqc

def load_config(path):
    with open(path, "r") as f:
        return yaml.safe_load(f)

def main():
    config = load_config("config/params.yaml")
    fastqc_cfg = config["fastqc"]

    dir_in = fastqc_cfg["dir_in"]
    dir_out = fastqc_cfg["dir_out"]
    dir_log = fastqc_cfg["dir_log"]
    file_format = fastqc_cfg["file_format"]
    thread = fastqc_cfg.get("thread", 4)
    job_name = fastqc_cfg["job_name"]
    job_limit = fastqc_cfg["job_limit"]
    partition = fastqc_cfg["partition"]

    run_fastqc(
        dir_in=dir_in,
        dir_out=dir_out,
        dir_log=dir_log,
        file_format=file_format,
        thread=thread,
        job_name=job_name,
        job_limit=job_limit,
        partition=partition,
    )

if __name__ == "__main__":
    main()