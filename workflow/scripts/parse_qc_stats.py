import re

def parse_qc_regex(rule, stat):
    ## Extract QC stats from logfile
    if rule == "whitelist":
        if stat=="n_reads_total":
            ## whitelist: Total reads
            return r"(\d+) (?:reads matched the barcode pattern)"
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

## featureCounts - summary - Assigned: Unique mapped, feature assigned, valid barcode & UMI reads.
##                                     equivalent to `n_reads_unimapped_valid - n_reads_skipped`
##                                     the denominator "n_reads" in the equation for calculating sequencing saturation
## STAR: Unique mapped, valid barcode & UMI reads

def get_stats(rule, stat, text):
    pattern = parse_qc_regex(rule, stat)
    m = re.search(pattern, text)
    return m.group(1)

logfile = "/home/luolab/GITHUB_REPO/ScRNAseq_smkpipe_at_Luolab/workflow/data/WangRuiyu/Example/logs/Testdata1/Testdata1_whitelist.log"
logfile = "/home/luolab/GITHUB_REPO/ScRNAseq_smkpipe_at_Luolab/workflow/data/YanTing/Pituitary/logs/YT032203/YT032203_whitelist.log"

rule = "whitelist"
stat = "n_reads_total"


## Two strategies:
##   1. re.findall + read(): speed varies, slower on large files
##   2. for loop + re.search: faster on large files
t=time.time();re.findall(parse_qc_regex("whitelist", "n_reads_total"), open(logfile).read());print(time.time()-t)
t=time.time();iter_log();print(time.time()-t)

re.findall(parse_qc_regex("whitelist", "n_reads_total"), open(logfile).read())




def iter_log():
    with open(logfile) as f:
        for line in f:
            pattern = parse_qc_regex(rule, stat)
            m = re.search(pattern, line)
            if m is not None:
                return m.group(1)