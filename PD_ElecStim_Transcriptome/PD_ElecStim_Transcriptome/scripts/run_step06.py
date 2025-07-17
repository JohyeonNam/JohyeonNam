### scripts/run_step06.py

##### Code by J. Nam
##### Last updated: 2025-07-17
##### Description: Run edgeR-based DE analysis from config

import subprocess
import yaml
from src.utils.mods_jnam import gtf_parser
import os

def load_config(path="config/params.yaml"):
    with open(path, "r") as f:
        return yaml.safe_load(f)

def main():
    config = load_config()
    edger_cfg = config["edger"]

    gtf_type = edger_cfg["gtf_type"]         # e.g., "GENCODE" or "NCBI"
    gtf_file = edger_cfg["gtf_path"]         # full path to GTF file
    geneinfo_outdir = os.path.dirname(edger_cfg["geneinfo_path"]) + "/"

    if not os.path.exists(os.path.join(geneinfo_outdir, f"custom_geneinfo_{gtf_type}.csv")):
        print("Generating geneinfo CSV from GTF...")
        gtf_parser(gtf_type, gtf_file, geneinfo_outdir)
    else:
        print(f"[INFO] Geneinfo file already exists. skipping generation.")

    cmd = [
        "Rscript", "src/step06_edger.R",
        edger_cfg["dir_in"],
        os.path.join(geneinfo_outdir, f"custom_geneinfo_{gtf_type}.csv"),
        edger_cfg["sampleinfo_path"],
        edger_cfg["dir_out"]
    ]

    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)

if __name__ == "__main__":
    main()