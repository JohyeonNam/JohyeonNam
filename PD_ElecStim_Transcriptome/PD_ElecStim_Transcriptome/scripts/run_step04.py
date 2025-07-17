#!/usr/bin/env python3

##### Code by J. Nam
##### Last updated: 2025-07-17
##### Description: Run STAR 2-pass alignment using config parameters

import yaml
from src.step04_star_align import run_star_alignment

def load_config(path):
    with open(path, "r") as f:
        return yaml.safe_load(f)

def main():
    config = load_config("config/params.yaml")
    align_cfg = config["star_align"]

    run_star_alignment(
        dir_in=align_cfg["dir_in"],
        dir_out=align_cfg["dir_out"],
        dir_log=align_cfg["dir_log"],
        file_format=align_cfg["file_format"],
        paired_end=align_cfg["paired_end"],
        thread=align_cfg.get("thread", 4),
        job_limit=align_cfg["job_limit"],
        partition=align_cfg["partition"],
        genome_index=align_cfg["genome_index"],
        job_name=align_cfg.get("job_name", "STARali")
    )

if __name__ == "__main__":
    main()