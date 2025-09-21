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
    visualise_prediction_cfg = config["visualise_prediction"]

    dir_out = visualise_prediction_cfg["dir_out"]
    dir_log = visualise_prediction_cfg["dir_log"]
    thread = visualise_prediction_cfg["thread"]
    job_limit = visualise_prediction_cfg["job_limit"]
    partition = visualise_prediction_cfg["partition"]
    job_name = visualise_prediction_cfg["job_name"]
    pred_df = visualise_prediction_cfg["pred_df"]

    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

    wrap_cmd = (
        f"time Rscript {project_root}/src/step16_visualise_prediction.R "
        f"{dir_out} {pred_df}"
    )
    
    subprocess.run(shlex.split(wrap_cmd), check=True)

if __name__ == "__main__":
    main()
