import os
import glob
import warnings
from pathlib import Path
import pandas as pd
from itertools import chain

# from snakemake.utils import validate

# validate(config, schema="../schemas/config.schema.yaml")

# ## Debug
# import yaml
# with open ('config/config_example.yaml') as f:
#     config = yaml.safe_load(f)

# Project Initialization
samples = (
    pd.read_csv(config["Samples"], dtype={'User': str, 'Project': str, 'Sample': str, 'File_R1': str, 'File_R2': str})
    .set_index("Sample", drop=False)
    .sort_index()
)

# samples.Samples must be unique
if samples["Sample"].duplicated().any():
    raise ValueError("Duplicated values found! The Sample column in samples_table.csv mustn't contain duplicated values.")

## Initialize project structure
Path('workflow', 'data', config['User'], config['Project']).mkdir(parents=True, exist_ok=True)

### Sub-dirs
sub_dirs = ['alignments', 'logs', 'miscs', 'outs']
sub_dir_config_keys = {'alignments':'aln_dir', 'logs':'log_dir', 'miscs':'misc_dir', 'outs':'out_dir'}

def init_paths(sub_dir):
    if not Path('workflow', 'data', config['User'], config['Project'], sub_dir).exists():
        if config[sub_dir_config_keys[sub_dir]] is None:
            print('\''+sub_dir_config_keys[sub_dir]+'\''+' not specified in config file. Creating \''+sub_dir+'\' in local folder.')
            Path('workflow', 'data', config['User'], config['Project'], sub_dir).mkdir(parents=True, exist_ok=True)
            config[sub_dir_config_keys[sub_dir]] = Path('workflow', 'data', config['User'], config['Project'], sub_dir)
        else:
            Path(config[sub_dir_config_keys[sub_dir]]).mkdir(parents=True, exist_ok=True)
            Path('workflow', 'data', config['User'], config['Project'], sub_dir).symlink_to(config[sub_dir_config_keys[sub_dir]])

for sub_dir in sub_dirs:
    init_paths(sub_dir)

### Source data dirs
if Path('workflow', 'data', config['User'], config['Project'], 'raw_fastqs').is_symlink():
    os.remove(str(Path('workflow', 'data', config['User'], config['Project'], 'raw_fastqs')))
if not Path('workflow', 'data', config['User'], config['Project'], 'raw_fastqs').exists():
    if config['src_fq_dir'] is None:
        raise ValueError('Missing input fastq files.')
    else:
        Path('workflow', 'data', config['User'], config['Project'], 'raw_fastqs').symlink_to(config['src_fq_dir'])
elif not Path('workflow', 'data', config['User'], config['Project'], 'raw_fastqs').is_symlink():
    config['src_fq_dir'] = Path('workflow', 'data', config['User'], config['Project'], 'raw_fastqs')

## 
def rget_fq_exts(src_dir, target, extensions):
    # fq = []
    # for ext in extensions:
        # fq = glob.glob(src_dir+'/**/'+target+'**'+ext, recursive=True)
        # fq.extend(Path(src_dir).rglob(target + "*" +ext))## pathlib.Path.rglob doesn't work for symlink
    fq = list(set(chain(*[glob.glob(src_dir+'/**/'+target+'/**'+ext, recursive=True) for ext in extensions])))
    if len(fq)>1:
        raise ValueError('More than one fastq found for target: ' + target + '! Execution halted.')
    elif len(fq)==0:
        raise ValueError('No fastq found for target: ' + target + '! Execution halted.')
    return str(fq[0])

def fetch_units_fq(units, config):
    fq_r1_extensions = config['r1_extensions']
    fq_r2_extensions = config['r2_extensions']
    for index, row in units.iterrows():
        src_dir = os.path.join("workflow", "data", row.User, row.Project, 'raw_fastqs')
        units.at[index,'path_to_R1'] = rget_fq_exts(src_dir, row['Library'], fq_r1_extensions)
        units.at[index,'path_to_R2'] = rget_fq_exts(src_dir, row['Library'], fq_r2_extensions)
        units.at[index,'File_R1'] = os.path.basename(units.at[index,'path_to_R1'])
        units.at[index,'File_R2'] = os.path.basename(units.at[index,'path_to_R2'])

# def fetch_table_fqs(table, table_type, config):
#     # If table is units, target='Library'; else if table is samples, target='Sample'
#     target = 'Library' if table_type=='units' else 'Sample' if table_type=='samples' else None
#     fq_r1_extensions = config['r1_extensions']
#     fq_r2_extensions = config['r2_extensions']
#     src_base = 'raw_fastqs' if target=='Library' else 'fastqs' if target == 'Sample' else None
#     for index, row in table.iterrows():
#         src_dir = os.path.join("workflow", "data", row.User, row.Project, src_base)
#         table.at[index,'path_to_R1'] = rget_fq_exts(src_dir, row[target], fq_r1_extensions)
#         table.at[index,'path_to_R2'] = rget_fq_exts(src_dir, row[target], fq_r2_extensions)

if config['Units'] is not None:
    # Case 1: some sample gets re-sequenced to gain deeper depth
    units = (
        pd.read_csv(config['Units'], dtype={'User': str, 'Project': str, 'Sample': str, 'Library': str, 'File_R1': str, 'File_R2': str})
        .set_index("Library", drop=False)
        .sort_index()
    )

    # units.Library must be unique
    # if not units["Library"].is_unique:
    if units["Library"].duplicated().any():
        raise ValueError("Duplicated values found! The Library column in units_table.csv mustn't contain duplicated values.")

    ## Dir to sample-level fastqs (aggregated units)
    Path('workflow', 'data', config['User'], config['Project'], 'fastqs').mkdir(parents=True, exist_ok=True)

    ## Glob source R1 R2 fq files, complete units table
    fetch_units_fq(units, config)

    ## Allocate space for aggregation of fq's
    if units['Sample'].duplicated().any():
        if config['tmp_fq_dir'] is None:
            raise ValueError('\'tmp_fq_dir\' not specified in config.yaml. Please specify a location to store aggregated fastq units. i.e. ' +
            str(Path('workflow', 'data', config['User'], config['Project'], 'tmp')))
        else:
            Path(config['tmp_fq_dir']).mkdir(parents=True, exist_ok=True)

    ## Prepare for aggregation (#Step 0 rule aggr_fqs)
    for sample in set(units['Sample']):
        u_df = units[units['Sample']==sample]
        if Path('workflow', 'data', config['User'], config['Project'], 'fastqs', sample).exists():
            os.remove(str(Path('workflow', 'data', config['User'], config['Project'], 'fastqs', sample)))
        if len(u_df)==1:
            path_to_sample = sorted(Path(config['src_fq_dir']).resolve().rglob(u_df['Library'].tolist()[0]+'/'))
            if len(path_to_sample) > 1:
                raise ValueError('More than one library observed under sample:' + sample + '. Execution halted.')
            Path('workflow', 'data', config['User'], config['Project'], 'fastqs', sample).symlink_to(path_to_sample[0])
            samples.at[sample, 'File_R1'] = units[units['Sample']==sample]["File_R1"][0]
            samples.at[sample, 'File_R2'] = units[units['Sample']==sample]["File_R2"][0]
        else:
            aggr_dir = Path(config['tmp_fq_dir'], config['User'], config['Project'], sample)
            aggr_dir.mkdir(parents=True, exist_ok=True)
            Path('workflow', 'data', config['User'], config['Project'], 'fastqs', sample).symlink_to(aggr_dir)
            samples.at[sample, 'File_R1'] = '_'.join([sample,'R1.fq.gz'])
            samples.at[sample, 'File_R2'] = '_'.join([sample,'R2.fq.gz'])
else:
    # Case 2: one sample per library, no re-sequencing




# ## Deprecated
# ## Parse sample/unit table, auto-infer R1 R2 files if not specified
# def fetch_fqR1R2(table, table_type, config):
#     fq_r1_extensions = config['r1_extensions']
#     fq_r2_extensions = config['r2_extensions']
#     for index, row in table.iterrows():
#         if table_type == "units":
#             tmp_dir = os.path.join("workflow", "data", row.User, row.Project, "raw_fastqs")
#         elif table_type == "samples":
#             tmp_dir = os.path.join("workflow", "data", row.User, row.Project, "fastqs")
#         ## TODO: Check that globbed file number equals to 1
#         if row.File_R1 not in os.listdir(tmp_dir):
#             table.at[index,'File_R1'] = os.path.basename(list(chain(*[glob.glob(os.path.join(tmp_dir, fq_r1_ext), recursive=True) for fq_r1_ext in fq_r1_extensions]))[0])
#         if row.File_R2 not in os.listdir(tmp_dir):
#             table.at[index,'File_R2'] = os.path.basename(list(chain(*[glob.glob(os.path.join(tmp_dir, fq_r2_ext), recursive=True) for fq_r2_ext in fq_r2_extensions]))[0])
#     return table
# 
# def fetch_fqR1R2(table, table_type, config):
#     fq_r1_extensions = config['r1_extensions']
#     fq_r2_extensions = config['r2_extensions']
#     for index, row in table.iterrows():
#         if table_type == "units":
#             src_dir = os.path.join("workflow", "data", row.User, row.Project, "raw_fastqs")
#             table.at[index,'File_R1'] = rget_fq_exts(src_dir, row.Library, fq_r1_extensions)
#             table.at[index,'File_R2'] = rget_fq_exts(src_dir, row.Library, fq_r2_extensions)
#         elif table_type == "samples":
#             src_dir = os.path.join("workflow", "data", row.User, row.Project, "fastqs")
#             table.at[index,'path_to_R1'] = rget_fq_exts(src_dir, row.Sample, fq_r1_extensions)
#             table.at[index,'path_to_R2'] = rget_fq_exts(src_dir, row.Sample, fq_r2_extensions)
#             table.at[index,'File_R1'] = os.path.basename(table.at[index,'path_to_R1'])
#             table.at[index,'File_R2'] = os.path.basename(table.at[index,'path_to_R2'])
#     return table
#
# for sample in set(units['Sample']):
#     u_df = units[units['Sample']==sample]
#     if Path('workflow', 'data', config['User'], config['Project'], 'fastqs', sample).exists():
#         os.remove(str(Path('workflow', 'data', config['User'], config['Project'], 'fastqs', sample)))
#     if len(u_df)==1:
#         path_to_sample = sorted(Path(config['src_fq_dir'].resolve()).rglob(u_df['Library'].tolist()[0]+'/'))
#         if len(path_to_sample) > 1:
#             raise ValueError('More than one library observed under sample:' + sample + '. Execution halted.')
#         Path('workflow', 'data', config['User'], config['Project'], 'fastqs', sample).symlink_to(path_to_sample[0])
#     else:
#         aggr_dir = Path(config['tmp_fq_dir'], config['User'], config['Project'], sample)
#         aggr_dir.mkdir(parents=True, exist_ok=True)
#         if not Path.joinpath(aggr_dir,sample+'_R1.fq.gz').exists():
#             cat_r1_command = u_df.File_R1.values.tolist()
#             cat_r1_command.insert(0,'cat')
#             cat_r1_command.append('> ' + str(aggr_dir) + '/' + sample + '_R1.fq.gz')
#             os.system(' '.join(cat_r1_command))
#         if not Path.joinpath(aggr_dir,sample+'_R2.fq.gz').exists():
#             cat_r2_command = u_df.File_R2.values.tolist()
#             cat_r2_command.insert(0,'cat')
#             cat_r2_command.append('> ' + str(aggr_dir) + '/' + sample + '_R2.fq.gz')
#             os.system(' '.join(cat_r2_command))
#         Path('workflow', 'data', config['User'], config['Project'], 'fastqs', sample).symlink_to(aggr_dir)
#
# samples = fetch_fqR1R2(samples, "samples", config)

def parse_suffix(rule):
    ## Given a rule name, return the suffix of its corresponding output file
    ## rule should be renamed, consider it later.
    if rule == 'umi_tools_whitelist':
        return 'whitelist.txt'
    elif rule == 'wash_whitelist':
        return 'whitelist_washed.txt'
    elif rule == 'umi_tools_extract':
        return 'extracted.fq.gz'
    elif rule == 'STAR':
        return 'Aligned.sortedByCoord.out.bam'
    elif rule == 'featureCounts':
        return 'Aligned.sortedByCoord.out.bam.featureCounts.bam'
    elif rule == 'sambamba_sort':
        return 'assigned_sorted.bam'
    elif rule == 'umi_tools_count':
        return 'counts_raw.tsv.gz'
    elif rule == "append_sfx":
        return 'counts.tsv.gz'
    elif rule == "umi_tools_whitelist_report":
        return 'cell_barcode_counts.png'
    elif rule == "featureCounts_report":
        return 'gene_assigned.summary'
    elif rule == "STAR_report":
        return 'Log.final.out'
    elif rule == "log_whitelist":
        return 'whitelist.log'
    elif rule == "log_extract":
        return 'extract.log'
    elif rule == "log_count":
        return 'count.log'
    elif rule == "velocyto":
        return 'velocyto.loom'

def get_table(tab):
    return units if tab=="units" else sample if tab=="samples" else None

def get_files(rule):
    files = expand(
        '_'.join(["workflow/data/{user}/{project}/alignments/{sample}/{sample}", parse_suffix(rule)]),
        zip, user=samples.User.to_list(), project=samples.Project.to_list(), sample=samples.Sample.to_list()
        )
    return files

def get_logfiles(rule):
    files = expand(
        '_'.join(["workflow/data/{user}/{project}/logs/{sample}/{sample}", parse_suffix(rule)]),
        zip, user=samples.User.to_list(), project=samples.Project.to_list(), sample=samples.Sample.to_list()
        )
    return files

def parse_dynamic_output(rule):
    # all_out: list of all files supposed to be generated in a given rule
    # Parse a list of output files. If a file was already produced in past run, remove that file from list.
    # Used to prevent redundant jobs on later reruns.
    all_out = get_files(rule)
    out_files = []
    for f in all_out: 
        if not os.path.exists(f):
            out_files.append(f)
    return out_files

def get_fastqs(wc):
    # Get input fastq read pairs from `samples` table
    s = samples.loc[wc.sample, ["User","Project","Sample","File_R1","File_R2"]].dropna()
    return {"read1": f"workflow/data/{s.User}/{s.Project}/fastqs/{s.Sample}/{s.File_R2}", 
            "read2": f"workflow/data/{s.User}/{s.Project}/fastqs/{s.Sample}/{s.File_R1}"}

def parse_STAR_dummy(wc):
    output = os.path.join("workflow", "data", wc.user, wc.project, "alignments", wc.sample, wc.sample) + "_Aligned.sortedByCoord.out.bam"
    if os.path.exists(output):
        return ""
    else:
        return "tmp/STARload.done"

def parse_fc_dummy(wc):
    output = os.path.join("workflow", "data", wc.user, wc.project, "alignments", wc.sample, wc.sample) + "_gene_assigned"
    if os.path.exists(output):
        return ""
    else:
        return "tmp/STARunload.done"

def get_report_output():
    return list(set(expand("workflow/data/{user}/{project}/outs/{project}_{suffix}", 
    zip, user=samples.User.to_list(), project=samples.Project.to_list(), suffix=['stats.csv','reads_stats.pdf','reads_stats.html'])))

def get_aggr_output():
    if config['RNA_velocity']:
        return list(set(expand("workflow/data/{user}/{project}/outs/{project}_{suffix}", 
        zip, user=samples.User.to_list(), project=samples.Project.to_list(), suffix=['counts_all.tsv.gz','velocyto_all.loom'])))
    else:
        return list(set(expand("workflow/data/{user}/{project}/outs/{project}_counts_all.tsv.gz", 
        zip, user=samples.User.to_list(), project=samples.Project.to_list())))
