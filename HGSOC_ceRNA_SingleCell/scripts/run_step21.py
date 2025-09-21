#!/usr/bin/python3

##### Code by J. Nam
##### Last updated: 2025-07-17
##### Description: Visualisation of the ORA results using parameters from config file

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
    plot_ora_cfg = config["plot_ora"]

    dir_out = plot_ora_cfg["dir_out"]
    dir_log = plot_ora_cfg["dir_log"]
    thread = plot_ora_cfg["thread"]
    job_limit = plot_ora_cfg["job_limit"]
    partition = plot_ora_cfg["partition"]
    job_name = plot_ora_cfg["job_name"]
    go_df = plot_ora_cfg["go_df"]
    kegg_df = plot_ora_cfg["kegg_df"]
    reactome_df = plot_ora_cfg["reactome_df"]
    wiki_df = plot_ora_cfg["wiki_df"]
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

    wrap_cmd = (
        f"time Rscript {project_root}/src/step21_plot_ora.R "
        f"{dir_out} {go_df} {kegg_df} {reactome_df} {wiki_df}"
    )
    
    subprocess.run(shlex.split(wrap_cmd), check=True)

if __name__ == "__main__":
    main()
