#!/usr/bin/python3

##### Code by J. Nam
##### Last updated: 2025-07-20
##### Code definition: Running CellRanger mkref

import glob
import os
import subprocess
from src.utils.mods_jnam import job_count_sleep, current_time

def run_cellranger_mkref(fasta_in, gtf_in, dir_out, dir_log, thread, job_limit, partition, job_name):
    start_time = current_time()

    if not fasta_in:
        raise FileNotFoundError(f"No {fasta_in} file found")
    if not gtf_in:
        raise FileNotFoundError(f"No {gtf_in} file found")
    
    file_name = os.path.join(dir_out, f"CellRanger_Ref_{start_time}")
    log_path = os.path.join(dir_log, f"{start_time}_CRmkref_{file_name}.log")

    wrap_cmd = (
        f"time /home/NJH/tools/cellranger-8.0.1/cellranger mkref --genome {file_name} "
        f"--fasta {fasta_in} --genes {gtf_in} --nthreads {thread}"
    )

    sbatch_cmd = [
        "sbatch", "-J", job_name, "-D", dir_log, "--output", log_path, "-p", partition,
        f'--wrap="{wrap_cmd}"'
    ]

    sbatch_str = " ".join(sbatch_cmd)
    print(f"Submitting job: {sbatch_str}")

    try:
        job_count_sleep(sbatch_str, job_name, int(job_limit) - 1)
    except Exception as e:
        print(f"[ERROR] Failed to submit job: {e}")