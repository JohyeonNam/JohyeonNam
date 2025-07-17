##### Code by J. Nam
##### Last updated: 2025-07-17
##### Code definition: Custom modules by J. Nam

import time
import subprocess
import os
import re

def current_time():
    """Return current time string in format YYMMDD_HHMMSS."""
    return time.strftime("%y%m%d_%H%M%S", time.localtime(time.time()))

def job_count_sleep(sbatch_cmd: str, job_name: str, job_limit: int):
    """Submit job using sbatch and wait if job count exceeds the limit."""
    try:
        subprocess.Popen(sbatch_cmd, shell=True)
        time.sleep(5)
        while subprocess.getoutput(f"squeue -n {job_name}").count(job_name) > job_limit:
            time.sleep(15)
    except Exception as e:
        raise RuntimeError(f"Job submission or monitoring failed: {e}")

def gtf_parser(gtf_type: str, gtf_file: str, out_path: str):
    """Extract GeneID, GeneName, and GeneType from GTF file and save as CSV."""
    output_lines = ["GeneID,GeneName,GeneType"]
    try:
        with open(gtf_file, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if len(fields) < 9 or fields[2] != "gene":
                    continue
                attr_str = fields[8]

                if gtf_type.upper() == "GENCODE":
                    geneid_match = re.search(r'gene_id "([^"]+)"', attr_str)
                    genename_match = re.search(r'gene_name "([^"]+)"', attr_str)
                    genetype_match = re.search(r'gene_type "([^"]+)"', attr_str)
                elif gtf_type.upper() == "NCBI":
                    geneid_match = re.search(r'db_xref "GeneID:(\d+)"', attr_str)
                    genename_match = re.search(r'gene_id "([^"]+)"', attr_str)
                    genetype_match = re.search(r'gene_biotype "([^"]+)"', attr_str)
                elif gtf_type.upper() == "ENSEMBL":
                    geneid_match = re.search(r'gene_id "([^"]+)"', attr_str)
                    genename_match = re.search(r'gene_name "([^"]+)"', attr_str) or geneid_match
                    genetype_match = re.search(r'gene_biotype "([^"]+)"', attr_str)
                else:
                    raise ValueError(f"Unsupported GTF type: {gtf_type}")

                if not (geneid_match and genename_match and genetype_match):
                    raise RuntimeError(f"Missing attribute in line: {attr_str}")

                geneid = geneid_match.group(1).split(".")[0] if gtf_type.upper() == "GENCODE" else geneid_match.group(1)
                genename = genename_match.group(1)
                genetype = genetype_match.group(1)
                output_lines.append(f"{geneid},{genename},{genetype}")
    except Exception as e:
        raise RuntimeError(f"Failed to parse GTF file: {e}")

    output_file = os.path.join(out_path, f"custom_geneinfo_{gtf_type}.csv")
    with open(output_file, "w") as f_out:
        f_out.write("\n".join(output_lines))