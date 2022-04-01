# import os
from pathlib import Path

user = snakemake.config["User"]
project = snakemake.config["Project"]
units = snakemake.input[0]

if snakemake.config['Units'] is not None:
    s_df = units[units['Sample']==snakemake.wildcards.sample]
    aggr_dir = Path(config['tmp_fq_dir'], config['User'], config['Project'], snakemake.wildcards.sample)
    
    # Combine r1
    cat_r1_command = s_df.path_to_R1.values.tolist()
    cat_r1_command.insert(0,'cat')
    cat_r1_command.append('> ' + str(aggr_dir) + '/' + snakemake.wildcards.sample + '_R1.fq.gz')
    os.system(' '.join(cat_r1_command))
    # Combine r2
    cat_r2_command = s_df.path_to_R2.values.tolist()
    cat_r2_command.insert(0,'cat')
    cat_r2_command.append('> ' + str(aggr_dir) + '/' + snakemake.wildcards.sample + '_R2.fq.gz')
    os.system(' '.join(cat_r2_command))
