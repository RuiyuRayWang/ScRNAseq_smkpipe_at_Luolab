import re
import pandas as pd
# import time

## Two strategies:
##   1. re.findall + read() + "big regex pattern": slower but can parse all stats in one line.
##   2. for loop + re.search + parse patterns individually: faster but more stats requires repeated function call, less readable
## Test runtime:
## If large amounts of stats need to be parsed per file, 1 is probably faster; If only a few stats need to be parsed per file, 2 is probably faster.
## With n_stats < 10, 2 is always faster 1.

# ## Slower than iter_log, mark as deprecated.
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

def iter_log(logfile, rule, stat):
    with open(logfile) as f:
        for line in f:
            pattern = parse_qc_regex(rule, stat)
            m = re.search(pattern, line)
            if m is not None:  ## Return as soon as a match is found. Dump the rest of the file.
                return m.group(1)

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
        elif stat=="n_BCs_filtered":
            ## extract: Filtered cell barcode. Not correctable.
            return r"^.*(?:Filtered cell barcode. Not correctable:) (\d+)$"
        elif stat=="n_BCs_corrected":
            ## extract: False cell barcode. Error-corrected.
            return r"^.*(?:False cell barcode. Error-corrected:) (\d+)$"
    elif rule == "count":
        if stat=="n_reads_unimapped_valid":
            ## count: Unique mapped, valid barcode & UMI reads
            ## equivalent to STAR "Uniquely mapped reads number"
            return r"^.*(?:Input Reads:) (\d+)$"
        elif stat=="n_reads_skipped":
            ## count: Unique mapped, BUT unassigned, valid barcode & UMI reads
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

print(snakemake.input.whitelist_log)
print(snakemake.input.whitelist_log[0])
print(type(snakemake.input.whitelist_log))




# stats = 


# logfile = "/home/luolab/GITHUB_REPO/ScRNAseq_smkpipe_at_Luolab/workflow/data/WangRuiyu/Example/logs/Testdata1/Testdata1_whitelist.log"
# logfile = "/home/luolab/GITHUB_REPO/ScRNAseq_smkpipe_at_Luolab/workflow/data/WangRuiyu/Example/logs/Testdata1/Testdata1_extract.log"
# logfile = "/home/luolab/GITHUB_REPO/ScRNAseq_smkpipe_at_Luolab/workflow/data/WangRuiyu/Example/logs/Testdata1/Testdata1_count.log"
# logfile = "/home/luolab/GITHUB_REPO/ScRNAseq_smkpipe_at_Luolab/workflow/data/YanTing/Pituitary/logs/YT032203/YT032203_whitelist.log"

# rule = "whitelist"
# stat = "n_reads_total"


# t=time.time();iter_log(logfile, rule="whitelist", stat="n_reads_total");iter_log(logfile, rule="whitelist", stat="n_BCs_unique");iter_log(logfile, rule="whitelist", stat="n_reads_selected");iter_log(logfile, rule="whitelist", stat="n_reads_correctable");print(time.time()-t)



# rule = "whitelist"
# stat = "n_reads_total"



# t=time.time();re.findall(parse_qc_regex("whitelist", "n_reads_total"), open(logfile).read());print(time.time()-t)
# t=time.time();re.findall(pattern, open(logfile).read());print(time.time()-t)
# t=time.time();iter_log();print(time.time()-t)


