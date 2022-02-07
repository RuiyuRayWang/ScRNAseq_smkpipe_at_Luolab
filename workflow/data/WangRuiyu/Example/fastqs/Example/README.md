# Example run

## Example Data
For example run, download the data from dropbox link below:
https://www.dropbox.com/sh/nwcqiqph3grhquw/AABuXp33_8VNXi8o8sByKSdYa?dl=1

## MD5SUM Check
Typically you'd expect to find an MD5 checksum file in the text format (i.e. MD5.txt), located in each of your library folder.

After the download is finished, md5sum of the example data can be checked against the following records:

```
ff360155afd18471e0debaa6dc98fc90  Example_R1.fq.gz
124cfc116c2820c66d8036fba737b8dc  Example_R2.fq.gz
```

## Execute

Activate a shell prompt, change to main directory of the repo (i.e. `ScRNAseq_smkpipe_at_Luolab/`) and execute the following command to start an example run:

```
snakemake --cores 32 --use-conda --configfile config/config_example.yaml
```
