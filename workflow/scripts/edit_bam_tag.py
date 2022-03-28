import os
import pysam

def parse_bam_tag(path_to_input_bam, path_to_output_bam):
    infile = pysam.AlignmentFile(path_to_input_bam, "rb")
    outfile = pysam.AlignmentFile(path_to_output_bam, "wb", template=infile)
    iter = infile.fetch(until_eof=True)
    for read in iter:
        read.set_tag('CB', str.split(read.query_name, sep='_')[1], replace=False)
        read.set_tag('UB', str.split(read.query_name, sep='_')[2], replace=False)
        read.query_name = str.split(read.query_name, sep='_')[0]
        outfile.write(read)
    infile.close()
    outfile.close()
    return None

parse_bam_tag(
    path_to_input_bam=os.path.abspath(snakemake.input[0]), 
    path_to_output_bam=os.path.abspath(snakemake.output[0])
    )
