#!/usr/bin/python3

##### Code by J. Nam
##### Last updated: 2025-07-17
##### Code definition: Running FastQC

import glob
import os
import subprocess
from src.utils.mods_jnam import job_count_sleep, current_time

def run_fastqc(dir_in, dir_out, dir_log, file_format, thread, job_limit, partition, job_name):
    start_time = current_time()
    fastq_list = glob.glob(os.path.join(dir_in, f"*.{file_format}"))
    if not fastq_list:
        dir_sublist = glob.glob(os.path.join(dir_in, "*"))
        for sub_dir in dir_sublist:
            fastq_list += glob.glob(os.path.join(sub_dir, f"*.{file_format}"))

    if not fastq_list:
        raise FileNotFoundError(f"No .{file_format} files found in {dir_in} or its subdirectories")

    total_jobs = len(fastq_list)
    current_job = 1

    for fastq_file in fastq_list:
        file_name = os.path.basename(fastq_file)
        log_path = os.path.join(dir_log, f"{start_time}_fastqc_{file_name}.log")

        wrap_cmd = (
            f'time /home/NJH/tools/FastQC/fastqc "{fastq_file}" '
            f'--extract -f fastq -t {thread} -o "{dir_out}"'
        )

        sbatch_cmd = [
            "sbatch", "-J", job_name, "-D", dir_log, "--output", log_path,
            "-p", partition, f'--wrap="{wrap_cmd}"'
        ]

        sbatch_str = " ".join(sbatch_cmd)
        print(f"Submitting job {current_job}/{total_jobs}: {sbatch_str}")

        try:
            job_count_sleep(sbatch_str, job_name, int(job_limit) - 1)
        except Exception as e:
            print(f"[ERROR] Failed to submit job: {e}")

        current_job += 1