#!/usr/bin/python3

##### Code by J. Nam
##### Last updated: 2025-07-17
##### Code definition: Running trimgalore

import glob
import os
from src.utils.mods_jnam import job_count_sleep, current_time

def run_trimgalore(dir_in, dir_out, dir_log, file_format, paired_end, thread, job_limit, partition, quality, job_name):
    start_time = current_time()

    if paired_end.upper() == 'P':
        fastq_list = glob.glob(os.path.join(dir_in, f"*1.{file_format}"))
        if not fastq_list:
            dir_sublist = glob.glob(os.path.join(dir_in, "*"))
            for sub_dir in dir_sublist:
                fastq_list += glob.glob(os.path.join(sub_dir, f"*1.{file_format}"))
    else:
        fastq_list = glob.glob(os.path.join(dir_in, f"*.{file_format}"))
        if not fastq_list:
            dir_sublist = glob.glob(os.path.join(dir_in, "*"))
            for sub_dir in dir_sublist:
                fastq_list += glob.glob(os.path.join(sub_dir, f"*.{file_format}"))

    if not fastq_list:
        raise FileNotFoundError(f"No matching .{file_format} files found in {dir_in} or its subdirectories")

    total_jobs = len(fastq_list)
    current_job = 1

    for fastq_file in fastq_list:
        file_name = os.path.basename(fastq_file)
        log_path = os.path.join(dir_log, f"{start_time}_trimgalore_{file_name}.log")

        if paired_end.upper() == 'P':
            reverse_file = fastq_file.replace(f"1.{file_format}", f"2.{file_format}")
            wrap_cmd = (
                f'time /home/NJH/tools/TrimGalore-0.6.7/trim_galore '
                f'--path_to_cutadapt /home/NJH/.local/bin/cutadapt '
                f'--paired -o "{dir_out}" --quality {quality} "{fastq_file}" "{reverse_file}"'
            )
        else:
            wrap_cmd = (
                f'time /home/NJH/tools/TrimGalore-0.6.7/trim_galore '
                f'--path_to_cutadapt /home/NJH/.local/bin/cutadapt '
                f'-o "{dir_out}" --quality {quality} "{fastq_file}"'
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
