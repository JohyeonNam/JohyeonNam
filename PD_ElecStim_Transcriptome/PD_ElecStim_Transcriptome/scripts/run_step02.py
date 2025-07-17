#!/usr/bin/env python3

##### Code by J. Nam
##### Last updated: 2025-07-17
##### Description: Run TrimGalore using parameters from config file

import yaml
from src.step02_trimgalore import run_trimgalore

def load_config(path):
    with open(path, "r") as f:
        return yaml.safe_load(f)

def main():
    config = load_config("config/params.yaml")
    trim_cfg = config["trimgalore"]

    run_trimgalore(
        dir_in=trim_cfg["dir_in"],
        dir_out=trim_cfg["dir_out"],
        dir_log=trim_cfg["dir_log"],
        file_format=trim_cfg["file_format"],
        paired_end=trim_cfg["paired_end"],
        thread=trim_cfg.get("thread", 4),
        job_limit=trim_cfg["job_limit"],
        partition=trim_cfg["partition"],
        quality=trim_cfg["quality"],
        job_name=trim_cfg.get("job_name", "TrimGal")
    )

if __name__ == "__main__":
    main()