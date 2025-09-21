#!/usr/bin/python3

##### Code by J. Nam
##### Last updated: 2025-07-20
##### Code definition: Running CellRanger count

import glob
import os
import subprocess
from src.utils.mods_jnam import job_count_sleep, current_time

def process_sample_sheet(sample_sheet):
    with open(sample_sheet, "r") as f:
        sample_info = f.readlines()
    
def run_cellranger_count(dir_in, dir_out, dir_log, thread, job_limit, partition, job_name, reference, sample_sheet):
    start_time = current_time()

    if not reference:
        raise FileNotFoundError(f"No {reference} file found")
    if not sample_sheet:
        raise FileNotFoundError(f"No {sample_sheet} file found")

    sample_info = process_sample_sheet(sample_sheet)
    
    total_jobs = len()
    current_job = 1

    for line in sample_info[1:]:
        try:
            pat_id = line.split(",")[0]
            runs = line.split(",")[2].split("/")
            runs = list(map(lambda x: dir_in+x, runs))
            fastq_files = ",".join(runs)

            log_path = os.path.join(dir_log, f"{start_time}_CRcount_{pat_id}.log")

            wrap_cmd = (
                f"time /home/NJH/tools/cellranger-8.0.1/cellranger count "
                f"--create-bam true --id={pat_id} --fastqs={fastq_files} "
                f"--transcriptome={reference} --sample={pat_id}"
            )

            sbatch_cmd = [
                "sbatch", "-J", job_name, "-D", dir_log, "--output", log_path, 
                "-p", partition, f'--wrap="{wrap_cmd}"'
            ]

            sbatch_str = " ".join(sbatch_cmd)
            print(f"Submitting CellRanger Count job {current_job}/{total_jobs}: {sbatch_str}")
            job_count_sleep(sbatch_str, job_name, int(job_limit) - 1)
            current_job += 1
    
    except Exception as e:
        print(f"[WARNING] Skipped {line}: {e}")