#!/usr/bin/python3

##### Code by J. Nam
##### Last updated: 2025-07-17
##### Code definition: Running STAR Align

import glob
import os
from src.utils.mods_jnam import job_count_sleep, current_time

def run_star_alignment(dir_in, dir_out, dir_log, file_format, paired_end, thread, job_limit, partition, genome_index, job_name):
    start_time = current_time()

    if paired_end.upper() == 'P':
        fastq_list = glob.glob(os.path.join(dir_in, f"*1_val_1.{file_format}"))
        if not fastq_list:
            raise FileNotFoundError("No paired-end FASTQ files found.")

        total_jobs = len(fastq_list)
        current_job = 1

        for fastq_file in fastq_list:
            reverse_file = fastq_file.replace(f"1_val_1.{file_format}", f"2_val_2.{file_format}")
            sample_name = os.path.basename(fastq_file).replace(f"1_val_1.{file_format}", "").rstrip("_")
            log_path = os.path.join(dir_log, f"{start_time}_STARalign_{sample_name}.log")
            prefix_path = os.path.join(dir_out, sample_name)

            wrap_cmd = (
                f"time /home/NJH/tools/STAR-2.7.10a/bin/Linux_x86_64/STAR "
                f"--runMode alignReads --outFilterType BySJout --limitBAMsortRAM 20000000000 "
                f"--outFilterMismatchNmax 999 --outFilterMultimapNmax 10 --alignSJoverhangMin 8 "
                f"--alignSJDBoverhangMin 1 --alignIntronMin 20 --alignIntronMax 1000000 "
                f"--alignMatesGapMax 1000000 --outFilterMismatchNoverLmax 0.02 "
                f"--runThreadN {thread} --genomeDir {genome_index} --readFilesIn {fastq_file} {reverse_file} "
                f"--outSAMtype BAM SortedByCoordinate --readFilesCommand zcat "
                f"--outWigType wiggle --outWigStrand Stranded --outWigNorm RPM "
                f"--twopassMode Basic --outFileNamePrefix {prefix_path}/"
            )

            sbatch_cmd = [
                "sbatch", "-J", job_name, "-D", dir_log, "--output", log_path,
                "-p", partition, f'--wrap="{wrap_cmd}"'
            ]

            sbatch_str = " ".join(sbatch_cmd)
            print(f"Submitting STAR alignment job {current_job}/{total_jobs}: {sbatch_str}")
            job_count_sleep(sbatch_str, job_name, int(job_limit) - 1)
            current_job += 1

    else:
        fastq_list = glob.glob(os.path.join(dir_in, f"*.{file_format}"))
        if not fastq_list:
            raise FileNotFoundError("No single-end FASTQ files found.")

        total_jobs = len(fastq_list)
        current_job = 1

        for fastq_file in fastq_list:
            sample_name = os.path.basename(fastq_file).replace(f".{file_format}", "").rstrip("_")
            log_path = os.path.join(dir_log, f"{start_time}_STARalign_{sample_name}.log")
            prefix_path = os.path.join(dir_out, sample_name)

            wrap_cmd = (
                f"time /home/NJH/tools/STAR-2.7.10a/bin/Linux_x86_64/STAR "
                f"--runMode alignReads --outFilterType BySJout --limitBAMsortRAM 20000000000 "
                f"--outFilterMismatchNmax 999 --outFilterMultimapNmax 10 --alignSJoverhangMin 8 "
                f"--alignSJDBoverhangMin 1 --alignIntronMin 20 --alignIntronMax 1000000 "
                f"--alignMatesGapMax 1000000 --outFilterMismatchNoverLmax 0.02 "
                f"--runThreadN {thread} --genomeDir {genome_index} --readFilesIn {fastq_file} "
                f"--outSAMtype BAM SortedByCoordinate --readFilesCommand zcat "
                f"--outWigType wiggle --outWigStrand Stranded --outWigNorm RPM "
                f"--twopassMode Basic --outFileNamePrefix {prefix_path}/"
            )

            sbatch_cmd = [
                "sbatch", "-J", job_name, "-D", dir_log, "--output", log_path,
                "-p", partition, f'--wrap="{wrap_cmd}"'
            ]

            sbatch_str = " ".join(sbatch_cmd)
            print(f"Submitting STAR alignment job {current_job}/{total_jobs}: {sbatch_str}")
            job_count_sleep(sbatch_str, job_name, int(job_limit) - 1)
            current_job += 1
