#!/usr/bin/python3

##### Code by J. Nam
##### Last updated: 2025-07-17
##### Code definition: Running featureCounts

import glob
import os
from src.utils.mods_jnam import job_count_sleep, current_time

def run_featurecounts(dir_in, dir_out, dir_log, paired_end, thread, job_limit, partition, job_name, feature_type, strand, gtf_path):
    start_time = current_time()

    paired_flag = " -p" if paired_end.upper() == "P" else ""
    sample_dirs = glob.glob(os.path.join(dir_in, "*"))

    total_jobs = len(sample_dirs)
    current_job = 1

    for sample_dir in sample_dirs:
        try:
            bam_files = glob.glob(os.path.join(sample_dir, "*.bam"))
            if not bam_files:
                continue

            bam_file = bam_files[0]
            sample_name = os.path.basename(sample_dir).replace("-", "_")
            log_path = os.path.join(dir_log, f"{start_time}_featureCounts_{sample_name}.log")
            output_prefix = os.path.join(dir_out, sample_name[:-8])

            wrap_cmd = (
                f"time /home/NJH/tools/subread-2.0.3-Linux-x86_64/bin/featureCounts "
                f"--countReadPairs -s {strand} -T {thread} -a {gtf_path} -t {feature_type}"
                f"{paired_flag} -o {output_prefix}.txt {bam_file}"
            )

            sbatch_cmd = [
                "sbatch", "-J", job_name, "-D", dir_log, "--output", log_path,
                "-p", partition, f'--wrap="{wrap_cmd}"'
            ]

            sbatch_str = " ".join(sbatch_cmd)
            print(f"Submitting featureCounts job {current_job}/{total_jobs}: {sbatch_str}")
            job_count_sleep(sbatch_str, job_name, int(job_limit) - 1)
            current_job += 1

        except Exception as e:
            print(f"[WARNING] Skipped {sample_dir}: {e}")