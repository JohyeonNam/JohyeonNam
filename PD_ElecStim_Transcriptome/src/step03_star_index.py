#!/usr/bin/python3

##### Code by J. Nam
##### Last updated: 2025-07-17
##### Code definition: Running STAR Index

import os
from src.utils.mods_jnam import job_count_sleep, current_time

def run_star_index(dir_log, fasta_path, gtf_path, thread, job_limit, partition, job_name):
    start_time = current_time()

    fasta_basename = os.path.basename(fasta_path)
    genome_dir = f"{fasta_basename}.index"
    log_path = os.path.join(dir_log, f"{start_time}_STARindex_{fasta_basename}.log")

    wrap_cmd = (
        f"time /home/NJH/tools/STAR-2.7.10a/bin/Linux_x86_64/STAR "
        f"--runMode genomeGenerate --genomeDir {genome_dir} "
        f"--genomeFastaFiles {fasta_path} --sjdbGTFfile {gtf_path} "
        f"--sjdbOverhang 100 --sjdbGTFfeatureExon exon --genomeSAindexNbases 14 "
        f"--runThreadN {thread}"
    )

    sbatch_cmd = [
        "sbatch", "-J", job_name, "-D", dir_log, "--output", log_path,
        "-p", partition, f'--wrap="{wrap_cmd}"'
    ]

    sbatch_str = " ".join(sbatch_cmd)
    print(f"Submitting STAR index job: {sbatch_str}")

    try:
        job_count_sleep(sbatch_str, job_name, int(job_limit) - 1)
    except Exception as e:
        print(f"[ERROR] Failed to submit STAR index job: {e}")
