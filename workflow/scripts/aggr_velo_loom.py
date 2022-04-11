import loompy
loompy.combine(snakemake.input, snakemake.output[0], key="Accession")