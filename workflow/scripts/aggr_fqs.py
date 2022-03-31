# import os
from pathlib import Path

user = snakemake.config["User"]
project = snakemake.config["Project"]

units = snakemake.input[0]

for sample in set(units['Sample']):
    u_df = units[units['Sample']==sample]
    aggr_dir = Path(config['tmp_fq_dir'], config['User'], config['Project'], sample)
    if len(u_df)>1:
        # Combine r1
        cat_r1_command = u_df.File_R1.values.tolist()
        cat_r1_command.insert(0,'cat')
        cat_r1_command.append('> ' + str(aggr_dir) + '/' + sample + '_R1.fq.gz')
        os.system(' '.join(cat_r1_command))
        # Combine r2
        cat_r2_command = u_df.File_R2.values.tolist()
        cat_r2_command.insert(0,'cat')
        cat_r2_command.append('> ' + str(aggr_dir) + '/' + sample + '_R2.fq.gz')
        os.system(' '.join(cat_r2_command))


