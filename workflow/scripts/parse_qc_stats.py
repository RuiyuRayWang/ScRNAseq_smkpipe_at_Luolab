import os
import re
import numpy as np
import pandas as pd
import altair as alt
# import seaborn as sns
# import matplotlib.pyplot as plt
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
        # if stat=="n_reads_total":  ## BUG: maxed-out at 100000001, producing negative values in n_reads_washed. Use extract.log to grab total input reads
        #     ## whitelist: Total reads
        #     return r"^.* (\d+) (?:reads matched the barcode pattern)$"
        if stat=="n_BCs_unique":
            ## whitelist: unique cell barcodes
            return r"^.*(?:Found) (\d+) (?:unique cell barcodes)$"
        elif stat=="n_reads_selected":
            ## whitelist: how many reads have a cell barcode in the whitelist
            return r"^.*(?:Found) (\d+) (?:total reads matching the selected cell barcodes)$"
        elif stat=="n_reads_correctable":
            ## whitelist: error correctable reads
            return r"^.*(?:Found) (\d+) (?:total reads which can be error corrected to the selected cell barcodes)$"
    elif rule == "extract":
        if stat=="n_reads_total":
            ## extract: Total reads
            return r"^.* (?:Input Reads:) (\d+)$"
        elif stat=="n_reads_valid_BC_UMI":
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
        elif stat=="n_reads_unassigned":
            ## count: Unique mapped, BUT unassigned, valid barcode & UMI reads
            ## Contains reads mapping to non-coding regions or intronic regions
            return r"^.*(?:Read skipped, no tag:) (\d+)$"
        elif stat=="n_reads_deduped":
            ## count: Unique mapped, feature assigned, valid barcode & UMI, deduplicated reads
            return r"^.*(?:Number of \(post deduplication\) reads counted:) (\d+)$"
## Note on more stats available but not parsed from the log:
## featureCounts - summary - Assigned: Unique mapped, feature assigned, valid barcode & UMI reads.
##                                     Name it "Effective Reads", denoted by n_reads_effective
##                                     n_reads_effective = n_reads_unimapped_valid - n_reads_unassigned`
##                                     Effective Reads is the denominator "n_reads" in the equation for calculating sequencing saturation.
## STAR: Unique mapped, valid barcode & UMI reads

def parse_logfiles(logfile, rule):
    l_name = os.path.basename(logfile)
    library = re.search("(.*)(?:_whitelist\.log|_count\.log|_extract\.log)", l_name).group(1)
    if rule == 'whitelist':
        return {
            'Library': library,
            'n_BCs_unique': iter_log_getstat(logfile, rule, 'n_BCs_unique'),
            'n_reads_selected': iter_log_getstat(logfile, rule, 'n_reads_selected'),
            'n_reads_correctable': iter_log_getstat(logfile, rule, 'n_reads_correctable')
        }
    elif rule == 'extract':
        return {
            'Library': library,
            'n_reads_total': iter_log_getstat(logfile, rule, 'n_reads_total'),
            'n_reads_valid_BC_UMI': iter_log_getstat(logfile, rule, 'n_reads_valid_BC_UMI'),
            'n_reads_uncorrectable': iter_log_getstat(logfile, rule, 'n_reads_uncorrectable'),
            'n_reads_corrected': iter_log_getstat(logfile, rule, 'n_reads_corrected')
        }
    elif rule == 'count':
        return {
            'Library': library,
            'n_reads_unimapped_valid': iter_log_getstat(logfile, rule, 'n_reads_unimapped_valid'),
            'n_reads_unassigned': iter_log_getstat(logfile, rule, 'n_reads_unassigned'),
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

## Parse stats from dictionary to dataframe
df_whitelist = pd.DataFrame(stats['whitelist'])
df_extract = pd.DataFrame(stats['extract'])
df_count = pd.DataFrame(stats['count'])

df_stats = pd.merge(pd.merge(df_whitelist, df_extract, how="left", on="Library"), df_count, how="left", on="Library").set_index('Library')
# df_stats["n_reads_washed"] = df_stats["n_reads_selected"] - df_stats["n_reads_valid_BC_UMI"]  ## BUG: output negative values, sth wrong, omitted
# df_stats["n_reads_filtered"] = df_stats["n_reads_uncorrectable"] - df_stats["n_reads_washed"]  ## Filtered becaused not selected
df_stats["n_reads_ununimapped"] = df_stats["n_reads_valid_BC_UMI"] - df_stats["n_reads_unimapped_valid"]  ## Multimapping + Unmapped too short + Unmapped other
df_stats["n_reads_effective"] = df_stats["n_reads_unimapped_valid"] - df_stats["n_reads_unassigned"]
df_stats["n_reads_collapsed"] = df_stats["n_reads_effective"] - df_stats["n_reads_deduped"]
df_stats["alignment_rate"] = df_stats["n_reads_unimapped_valid"] / df_stats["n_reads_valid_BC_UMI"]
df_stats["seq_saturation"] = 1 - df_stats["n_reads_deduped"] / df_stats["n_reads_effective"]
df_stats.to_csv(snakemake.output[0], index=True, sep=",", header=True)

# ## Plots
# ### Reads perspective
# df_plot = df_stats[["n_reads_deduped","n_reads_collapsed","n_reads_unassigned","n_reads_ununimapped","n_reads_uncorrectable"]]
# labs = ["Final Reads", "Deduplicated Reads", "Reads Not Assigned to Feature", "Not Unique Mapped Reads", "Non-correctable Reads"]
# colors = sns.color_palette("Set2", 5)

# fig, ax = plt.subplots()

# left = np.zeros(len(df_plot))

# for i, col in enumerate(df_plot.columns):
#     ax.barh(
#         df_plot.index, df_plot[col], left=left, label=labs[i], height=0.36, edgecolor="#000000", color=colors[i])
#     left += np.array(df_plot[col])

# totals = df_plot.sum(axis=1)
# for i, total in enumerate(totals):
#     ax.text(total, totals.index[i], round(total), va='center',
#             weight='bold', rotation=270)

# x_offset = - max(totals) * 0.03
# for bar in ax.patches:
#     ax.text(
#         bar.get_width() + bar.get_x() + x_offset,
#         bar.get_y() + bar.get_height() / 2,
#         round(bar.get_width()),
#         va='center',
#         color='w',
#         weight='bold',
#         size=8,
#         rotation=270
#     )

# ax.set_title('QC Statistics: Reads Perspective')
# ax.set_xlabel('# Reads')
# ax.set_ylabel('Library')

# box = ax.get_position()
# ax.set_position([box.x0, box.y0 + box.height * 0.1,
#                  box.width, box.height * 0.9])

# ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
#           fancybox=True, shadow=True, ncol=2)

# # fig.tight_layout()
# # fig.show()
# fig.savefig(snakemake.output[1], dpi=150, bbox_inches="tight")

## Bar chart with interactive legend
df_plot = df_stats[["n_reads_deduped","n_reads_collapsed","n_reads_unassigned","n_reads_ununimapped","n_reads_uncorrectable"]]
df_plot = df_plot.melt(var_name='stats', value_name='n_reads', ignore_index=False).reset_index()
df_plot['order'] = df_plot['stats'].replace(
    {val: i for i, val in enumerate(["n_reads_deduped","n_reads_collapsed","n_reads_unassigned","n_reads_ununimapped","n_reads_uncorrectable"])}
)

selection = alt.selection_multi(fields=['stats'], bind='legend')

bars = alt.Chart(df_plot).mark_bar().encode(
    x=alt.X('sum(n_reads):Q', stack='zero'),
    y=alt.Y('Library:N'),
    color=alt.Color('stats', 
        sort=alt.EncodingSortField('order', order='descending')
    ),
    order='order',
    opacity=alt.condition(selection, alt.value(1), alt.value(0.2))
).add_selection(
    selection
)

bars.save(snakemake.output[1])
bars.save(snakemake.output[2])

# =================== DEBUG ===================

# ## For debug purposes
# log_whitelist = ['workflow/data/WangRuiyu/Example/logs/Testdata1/Testdata1_whitelist.log', 'workflow/data/WangRuiyu/Example/logs/Testdata2/Testdata2_whitelist.log']
# log_extract = ['workflow/data/WangRuiyu/Example/logs/Testdata1/Testdata1_extract.log', 'workflow/data/WangRuiyu/Example/logs/Testdata2/Testdata2_extract.log']
# log_count = ['workflow/data/WangRuiyu/Example/logs/Testdata1/Testdata1_count.log', 'workflow/data/WangRuiyu/Example/logs/Testdata2/Testdata2_count.log']

# import glob

# log_whitelist = glob.glob("workflow/data/YanTing/Pituitary/logs/**/*whitelist.log", recursive=True)
# log_extract = glob.glob("workflow/data/YanTing/Pituitary/logs/**/*extract.log", recursive=True)
# log_count = glob.glob("workflow/data/YanTing/Pituitary/logs/**/*count.log", recursive=True)

# for rule in rules:
#     if rule == 'whitelist':
#         for lf in log_whitelist:
#             stats[rule].append(parse_logfiles(lf, rule))
#     elif rule == 'extract':
#         for lf in log_extract:
#             stats[rule].append(parse_logfiles(lf, rule))
#     elif rule == 'count':
#         for lf in log_count:
#             stats[rule].append(parse_logfiles(lf, rule))

# ================ DEPREACATED ================

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
#     "count":    ["n_reads_unimapped_valid","n_reads_unassigned","n_reads_deduped"]
#     }

## runtime test for strategy [1]
# t=time.time();iter_log_getstat(logfile, rule="whitelist", stat="n_reads_total");iter_log_getstat(logfile, rule="whitelist", stat="n_BCs_unique");iter_log_getstat(logfile, rule="whitelist", stat="n_reads_selected");iter_log_getstat(logfile, rule="whitelist", stat="n_reads_correctable");print(time.time()-t)

### Deprecated Plotting Mechanics
## Plot test
# ax = df_stats[["n_reads_deduped","n_reads_collapsed","n_reads_unassigned","n_reads_ununimapped","n_reads_uncorrectable"]].plot.barh(stacked=True)
# plt.show()

## Plot horizontal
# df_test = df_stats[["n_reads_deduped","n_reads_collapsed"]]

# fig, ax = plt.subplots()

# colors = ['#24b1d1', '#ae24d1']
# bottom = np.zeros(len(df_test))

# for i, col in enumerate(df_test.columns):
#     ax.bar(
#         df_test.index, df_test[col], bottom=bottom, label=col, color=colors[i])
#     bottom += np.array(df_test[col])

# totals = df_test.sum(axis=1)
# y_offset = 10000
# for i, total in enumerate(totals):
#     ax.text(totals.index[i], total + y_offset, round(total), ha='center',
#             weight='bold')

# y_offset = -25000
# for bar in ax.patches:
#     ax.text(
#         bar.get_x() + bar.get_width() / 2,
#         bar.get_height() + bar.get_y() + y_offset,
#         round(bar.get_height()),
#         ha='center',
#         color='w',
#         weight='bold',
#         size=8
#     )

# ax.set_title('Tips by Day and Gender')
# ax.legend()
# plt.show()

# fig.savefig("/home/luolab/GITHUB_REPO/test.png", dpi=150, bbox_inches="tight")