import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# import time

## Two strategies:
##   [1]. re.findall + read() + "big regex pattern": slower but can parse all stats in one line.
##   [2]. for loop + re.search + parse patterns individually: faster but more stats requires repeated function call, less readable
## Tested runtime:
## If large amounts of stats need to be parsed per file, [1] is probably faster; If only a few stats need to be parsed per file, [2] is probably faster.
## Runtime tests show that with n_stats < 10, [2] is always faster [1].

def iter_log_getstat(logfile, rule, stat):
    with open(logfile) as f:
        for line in f:
            pattern = parse_qc_regex(rule, stat)
            m = re.search(pattern, line)
            if m is not None:  ## Return as soon as a match is found. Dump the rest of the file.
                return int(m.group(1))

def parse_qc_regex(rule, stat):
    ## Extract QC stats from logfile
    if rule == "whitelist":
        if stat=="n_reads_total":
            ## whitelist: Total reads
            return r"^.* (\d+) (?:reads matched the barcode pattern)$"
        elif stat=="n_BCs_unique":
            ## whitelist: unique cell barcodes
            return r"^.*(?:Found) (\d+) (?:unique cell barcodes)$"
        elif stat=="n_reads_selected":
            ## whitelist: 
            return r"^.*(?:Found) (\d+) (?:total reads matching the selected cell barcodes)$"
        elif stat=="n_reads_correctable":
            ## whitelist: error correctable reads
            return r"^.*(?:Found) (\d+) (?:total reads which can be error corrected to the selected cell barcodes)$"
    elif rule == "extract":
        if stat=="n_reads_valid_BC_UMI":
            ## extract: reads with valid barcode & UMI
            return r"^.*(?:Reads output:) (\d+)$"
        elif stat=="n_reads_uncorrectable":
            ## extract: Filtered cell barcode. Not correctable.
            return r"^.*(?:Filtered cell barcode. Not correctable:) (\d+)$"
        elif stat=="n_reads_corrected":
            ## extract: False cell barcode. Error-corrected.
            return r"^.*(?:False cell barcode. Error-corrected:) (\d+)$"
    elif rule == "count":
        if stat=="n_reads_unimapped_valid":
            ## count: Unique mapped, valid barcode & UMI reads
            ## equivalent to STAR "Uniquely mapped reads number"
            return r"^.*(?:Input Reads:) (\d+)$"
        elif stat=="n_reads_skipped":
            ## count: Unique mapped, BUT unassigned, valid barcode & UMI reads
            ## Contains reads mapping to non-coding regions or intronic regions
            return r"^.*(?:Read skipped, no tag:) (\d+)$"
        elif stat=="n_reads_deduped":
            ## count: Unique mapped, feature assigned, valid barcode & UMI, deduplicated reads
            return r"^.*(?:Number of \(post deduplication\) reads counted:) (\d+)$"
## Note on more stats available but not parsed from the log:
## featureCounts - summary - Assigned: Unique mapped, feature assigned, valid barcode & UMI reads.
##                                     Name it "Effective Reads", denoted by n_reads_effective
##                                     n_reads_effective = n_reads_unimapped_valid - n_reads_skipped`
##                                     Effective Reads is the denominator "n_reads" in the equation for calculating sequencing saturation.
## STAR: Unique mapped, valid barcode & UMI reads

def parse_logfiles(logfile, rule):
    library = re.search("^.*/([^_]*)_([^_]*).log$", logfile).group(1)
    if rule == 'whitelist':
        return {
            'Library': library,
            'n_reads_total': iter_log_getstat(logfile, rule, 'n_reads_total'),
            'n_BCs_unique': iter_log_getstat(logfile, rule, 'n_BCs_unique'),
            'n_reads_selected': iter_log_getstat(logfile, rule, 'n_reads_selected'),
            'n_reads_correctable': iter_log_getstat(logfile, rule, 'n_reads_correctable')
        }
    elif rule == 'extract':
        return {
            'Library': library,
            'n_reads_valid_BC_UMI': iter_log_getstat(logfile, rule, 'n_reads_valid_BC_UMI'),
            'n_reads_uncorrectable': iter_log_getstat(logfile, rule, 'n_reads_uncorrectable'),
            'n_reads_corrected': iter_log_getstat(logfile, rule, 'n_reads_corrected')
        }
    elif rule == 'count':
        return {
            'Library': library,
            'n_reads_unimapped_valid': iter_log_getstat(logfile, rule, 'n_reads_unimapped_valid'),
            'n_reads_skipped': iter_log_getstat(logfile, rule, 'n_reads_skipped'),
            'n_reads_deduped': iter_log_getstat(logfile, rule, 'n_reads_deduped')
        }

## Define rules for which logs will be parsed
rules = ['whitelist', 'extract', 'count']

## Instantiate an empty dict to store stats
stats = {'whitelist':[], 'extract':[], 'count':[]}

for rule in rules:
    if rule == 'whitelist':
        for lf in snakemake.input.log_whitelist:
            stats[rule].append(parse_logfiles(lf, rule))
    elif rule == 'extract':
        for lf in snakemake.input.log_extract:
            stats[rule].append(parse_logfiles(lf, rule))
    elif rule == 'count':
        for lf in snakemake.input.log_count:
            stats[rule].append(parse_logfiles(lf, rule))

for rule in rules:
    if rule == 'whitelist':
        for lf in log_whitelist:
            stats[rule].append(parse_logfiles(lf, rule))
    elif rule == 'extract':
        for lf in log_extract:
            stats[rule].append(parse_logfiles(lf, rule))
    elif rule == 'count':
        for lf in log_count:
            stats[rule].append(parse_logfiles(lf, rule))            

## Parse stats from dictionary to dataframe
df_whitelist = pd.DataFrame(stats['whitelist'])
df_extract = pd.DataFrame(stats['extract'])
df_count = pd.DataFrame(stats['count'])

df_stats = pd.merge(pd.merge(df_whitelist, df_extract, how="left", on="Library"), df_count, how="left", on="Library").set_index('Library')
df_stats["n_reads_washed"] = df_stats["n_reads_selected"] - df_stats["n_reads_valid_BC_UMI"]
df_stats["n_reads_ununimapped"] = df_stats["n_reads_valid_BC_UMI"] - df_stats["n_reads_unimapped_valid"]  ## Multimapping + Unmapped too short + Unmapped other
df_stats["n_reads_effective"] = df_stats["n_reads_unimapped_valid"] - df_stats["n_reads_skipped"]
df_stats["n_reads_collapsed"] = df_stats["n_reads_effective"] - df_stats["n_reads_deduped"]
df_stats["alignment_rate"] = df_stats["n_reads_unimapped_valid"] / df_stats["n_reads_valid_BC_UMI"]
df_stats["seq_saturation"] = 1 - df_stats["n_reads_deduped"] / df_stats["n_reads_effective"]
df_stats.to_csv(snakemake.output[0], index=True, sep=",", header=True)


# ================ DEPREACATED ================

## For debug


# ## Strategy [2]
# ## Slower than iter_log_getstat, marked as deprecated.
# ## whitelist big pattern
# pattern = r".* (\d+) (?:reads matched the barcode pattern\n)"+\
#     r".*(?:Found) (\d+) (?:unique cell barcodes\n)"+\
#     r".*(?:Found) (\d+) (?:total reads matching the selected cell barcodes\n)"+\
#     r".*(?:Found) (\d+) (?:total reads which can be error corrected to the selected cell barcodes)"

# ## extract big pattern
# pattern = r".*(?:Reads output:) (\d+)\n"+\
#     r".*(?:Filtered cell barcode. Not correctable:) (\d+)\n"+\
#     r".*(?:False cell barcode. Error-corrected:) (\d+)"

# ## count big pattern
# pattern = r".*(?:Input Reads:) (\d+)\n"+\
#     r".*(?:Read skipped, no tag:) (\d+)\n"+\
#     r".*(?:Number of \(post deduplication\) reads counted:) (\d+)"

# t=time.time();re.findall(pattern, open(logfile).read());print(time.time()-t)
# pd.DataFrame(m, columns=['stat1','stat2','stat3']);print(time.time()-t)


## equality check
# len(snakemake.input.log_whitelist) == len(snakemake.input.log_extract) == len(snakemake.input.log_count)

# rule_stats = {
#     "whitelist":["n_reads_total","n_BCs_unique","n_reads_selected","n_reads_correctable"],
#     "extract":  ["n_reads_valid_BC_UMI","n_reads_uncorrectable","n_reads_corrected"],
#     "count":    ["n_reads_unimapped_valid","n_reads_skipped","n_reads_deduped"]
#     }


## runtime test for strategy [1]
# t=time.time();iter_log_getstat(logfile, rule="whitelist", stat="n_reads_total");iter_log_getstat(logfile, rule="whitelist", stat="n_BCs_unique");iter_log_getstat(logfile, rule="whitelist", stat="n_reads_selected");iter_log_getstat(logfile, rule="whitelist", stat="n_reads_correctable");print(time.time()-t)
