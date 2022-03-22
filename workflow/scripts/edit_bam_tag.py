"""
USAGE:-

python edit_bam_tag.py -i <input_file> -o <output_file>

"""

import os
import argparse
import pysam

parser = argparse.ArgumentParser(description = 'Generate BAM with OQ tag')
parser.add_argument('-i', '--input', required=True, help='Input mark dup BAM File')
parser.add_argument('-o', '--output', required=True, help='Output BAM file with OQ tag')

args = parser.parse_args()
infile_path = os.path.abspath(args.input)
outfile_path = os.path.abspath(args.output)

infile = pysam.AlignmentFile(infile_path, "rb")
outfile = pysam.AlignmentFile(outfile_path, "wb", template=infile)
iter = infile.fetch(until_eof=True)
for read in iter:
    read.set_tag('CB', str.split(read.query_name, sep='_')[1], replace=False)
    read.set_tag('UB', str.split(read.query_name, sep='_')[2], replace=False)
    read.query_name = str.split(read.query_name, sep='_')[0]
    outfile.write(read)
infile.close()
outfile.close()