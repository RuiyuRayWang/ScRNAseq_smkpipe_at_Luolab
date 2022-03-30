# import os
from pathlib import Path

user=snakemake.config["User"]
project=snakemake.config["Project"]

Path('workflow','data',user,project).mkdir(parents=True, exist_ok=True)

# alignments
if snakemake.config['aln_dir'] is None:
    Path('workflow','data',user,project,'alignment').mkdir(parents=True, exist_ok=True)
else:
    Path('workflow','data',user,project,'alignment').symlink_to(snakemake.config['aln_dir'])

# logs
if snakemake.config['log_dir'] is None:
    Path('workflow','data',user,project,'logs').mkdir(parents=True, exist_ok=True)
else:
    Path('workflow','data',user,project,'logs').symlink_to(snakemake.config['log_dir'])

# miscs
if snakemake.config['misc_dir'] is None:
    Path('workflow','data',user,project,'miscs').mkdir(parents=True, exist_ok=True)
else:
    Path('workflow','data',user,project,'miscs').symlink_to(snakemake.config['misc_dir'])

# outs
if snakemake.config['out_dir'] is None:
    Path('workflow','data',user,project,'outs').mkdir(parents=True, exist_ok=True)
else:
    Path('workflow','data',user,project,'outs').symlink_to(snakemake.config['out_dir'])

# fastqs
if snakemake.config['src_fq_dir'] is None:
    if not Path('workflow','data',user,project,'raw_fastqs').exists():
        raise ValueError("Missing Input fastq files.")
else:
    Path('workflow','data',user,project,'raw_fastqs').symlink_to(snakemake.config['src_fq_dir'])

# parse samples and libraries

