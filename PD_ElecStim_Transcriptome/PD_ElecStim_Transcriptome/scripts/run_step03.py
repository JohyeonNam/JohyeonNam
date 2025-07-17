#!/usr/bin/env python3

##### Code by J. Nam
##### Last updated: 2025-07-17
##### Description: Run STAR index generation using config parameters

import yaml
from src.step03_star_index import run_star_index

def load_config(path):
    with open(path, "r") as f:
        return yaml.safe_load(f)

def main():
    config = load_config("config/params.yaml")
    star_cfg = config["star_index"]

    run_star_index(
        dir_log=star_cfg["dir_log"],
        fasta_path=star_cfg["fasta_path"],
        gtf_path=star_cfg["gtf_path"],
        thread=star_cfg.get("thread", 4),
        job_limit=star_cfg["job_limit"],
        partition=star_cfg["partition"],
        job_name=star_cfg.get("job_name", "STARind")
    )

if __name__ == "__main__":
    main()