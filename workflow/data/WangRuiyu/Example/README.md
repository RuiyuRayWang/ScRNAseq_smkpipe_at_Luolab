# Example run

## Example Data
For example run, download the data from dropbox link below:
https://www.dropbox.com/sh/nwcqiqph3grhquw/AABuXp33_8VNXi8o8sByKSdYa?dl=1

We provide three libraries (down-sampled) from two batches as example datasets. Each of these libraries contains two fastq files in the format of gzip compressed `R1`, `R2` pairs.

## MD5SUM Check

After the download is finished, md5sum's of the example dataset can be checked against the following records.   

```
# Testdata1_1
54381921e4984307a0b742bb454b7a20  Testdata1_1_R1.fq.gz
ed9e05592f15af2a94c1e5208d57e207  Testdata1_1_R2.fq.gz

# Testdata2
5acaddaf543988ba9fe1b1afb112ca70  Testdata2_R1.fq.gz
3025f8a97b67ee5c8624e8d61095c198  Testdata2_R2.fq.gz

# Testdata1_2
88f1f10960adef5665af96247025ef43  Testdata1_2_R1.fq.gz
7da534bfecb5c129e946e12b81d50eef  Testdata1_2_R2.fq.gz
```

Alternatively, you may run a checksum against the `MD5.txt` files located at `raw_fastqs/Testdata*/MD5.txt`:

```
$ cd raw_fastqs/Batch_1/2.cleandata/Testdata1_1
$ md5sum -c MD5.txt > checkresult
```

## Execute

Activate a shell prompt, navigate to the main directory of the repository (i.e. `ScRNAseq_smkpipe_at_Luolab/`) and execute the following command to start an example run:

```
snakemake --cores 32 --use-conda --configfile config/config_example.yaml
```
