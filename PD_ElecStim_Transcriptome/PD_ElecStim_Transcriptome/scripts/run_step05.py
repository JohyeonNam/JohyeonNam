#!/usr/bin/env python3

##### Code by J. Nam
##### Last updated: 2025-07-17
##### Description: Run featureCounts using config parameters

import yaml
from src.step05_featurecounts import run_featurecounts

def load_config(path):
    with open(path, "r") as f:
        return yaml.safe_load(f)

def main():
    config = load_config("config/params.yaml")
    fc_cfg = config["featurecounts"]

    run_featurecounts(
        dir_in=fc_cfg["dir_in"],
        dir_out=fc_cfg["dir_out"],
        dir_log=fc_cfg["dir_log"],
        paired_end=fc_cfg["paired_end"],
        thread=fc_cfg.get("thread", 4),
        job_limit=fc_cfg["job_limit"],
        partition=fc_cfg["partition"],
        job_name=fc_cfg.get("job_name", "featCnt"),
        feature_type=fc_cfg["feature_type"],
        strand=fc_cfg["strand"],
        gtf_path=fc_cfg["gtf_path"]
    )

if __name__ == "__main__":
    main()