# Example run

## Example Data
For example run, download the data from dropbox link below:
https://www.dropbox.com/sh/nwcqiqph3grhquw/AABuXp33_8VNXi8o8sByKSdYa?dl=1

We provide two libraries as an example dataset, each contains two fastq files in the format of gzip compressed `R1`, `R2` pairs.

## MD5SUM Check
You should find an MD5 checksum file, in the text format (i.e. MD5.txt), located in each of the library folder.

After the download is finished, md5sum's of the example dataset can be checked against the following records.  
You can find the same md5sums in the `fastqs/Testdata*/MD5.txt` files in this repo:

```
# Testdata1
ff360155afd18471e0debaa6dc98fc90  Testdata1_R1.fq.gz
124cfc116c2820c66d8036fba737b8dc  Testdata1_R2.fq.gz

# Testdata2
5acaddaf543988ba9fe1b1afb112ca70  Testdata2_R1.fq.gz
3025f8a97b67ee5c8624e8d61095c198  Testdata2_R2.fq.gz
```

## Execute

Activate a shell prompt, navigate to the main directory of the repository (i.e. `ScRNAseq_smkpipe_at_Luolab/`) and execute the following command to start an example run:

```
snakemake --cores 32 --use-conda --configfile config/config_example.yaml
```
