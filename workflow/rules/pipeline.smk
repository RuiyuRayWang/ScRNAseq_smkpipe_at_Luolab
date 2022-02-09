## TODO: Initialize project structure per config.yaml - user/project/etc.

## TODO: FASTQ file quality control & summarize QC statistics

# Step 1: Identify cell barcode whitelist (identify correct BC)
rule umi_tools_whitelist:
    input:
        unpack(get_fastqs),
    output:
        "workflow/data/{user}/{project}/alignments/{library}/{library}_whitelist.txt"
    log:
        "workflow/data/{user}/{project}/logs/{library}/{library}_whitelist.log"
    conda:
        "../envs/master.yaml"
    threads:
        1
    shell:
        """
        umi_tools whitelist --bc-pattern=CCCCCCCCNNNNNNNN \
                            --log={log} \
                            --stdin {input.read1} \
                            --set-cell-number=100 \
                            --plot-prefix=workflow/data/{wildcards.user}/{wildcards.project}/alignments/{wildcards.library}/{wildcards.library} \
                            --log2stderr > {output}
        """

# Step 2: Wash whitelist
rule wash_whitelist:
    input:
        whitelist="workflow/data/{user}/{project}/alignments/{library}/{library}_whitelist.txt",
        ground_truth=config['ground_truth']
    output:
        "workflow/data/{user}/{project}/alignments/{library}/{library}_whitelist_washed.txt"
    threads:
        1
    script:
        "../scripts/wash_whitelist.py"

# Step 3: Extract barcodes and UMIs and add to read names
rule umi_tools_extract:
    input:
        unpack(get_fastqs),
        whitelist_washed="workflow/data/{user}/{project}/alignments/{library}/{library}_whitelist_washed.txt"
    output:
        "workflow/data/{user}/{project}/alignments/{library}/{library}_extracted.fq.gz"
    log:
        "workflow/data/{user}/{project}/logs/{library}/{library}_extract.log"
    conda:
        "../envs/master.yaml"
    threads:
        1
    shell:
        """
        umi_tools extract --bc-pattern=CCCCCCCCNNNNNNNN \
                          --log={log} \
                          --stdin {input.read1} \
                          --read2-in {input.read2} \
                          --stdout {output} \
                          --read2-stdout \
                          --filter-cell-barcode \
                          --error-correct-cell \
                          --whitelist={input.whitelist_washed}
        """

# # Step 4-0: Generate genome index
rule STAR_gen:
    input:
        genome_fa=config["genome_fa"],
        gtf=config["gtf_annotation"]
    output:
        directory(config["genome_index"])
    params:
        sjdbOverhang=config["sjdbOverhang"]
    conda:
        "../envs/master.yaml"
    threads:
        config["threads"]
    shell:
        """
        STAR --runThreadN {threads} \
             --runMode genomeGenerate \
             --genomeDir {output} \
             --genomeFastaFiles {input.genome_fa} \
             --sjdbGTFfile {input.gtf} \
             --sjdbOverhang {params.sjdbOverhang}
        """

# # Step 4-1: Load Genome Index
rule STAR_load:
    input:
        ex=get_files('umi_tools_extract'),
        genomeDir=config["genome_index"]
    output:
        temp(touch("tmp/STARload.done"))
    conda:
        "../envs/master.yaml"
    threads:
        config["threads"]
    shell:
        """
        STAR --runThreadN {threads} \
             --genomeLoad LoadAndExit \
             --genomeDir {input.genomeDir} \
             --outSAMmode None
        """

rule STAR:
    input:
        extracted_fq="workflow/data/{user}/{project}/alignments/{library}/{library}_extracted.fq.gz",
        genomeDir=config["genome_index"],
        dummy=parse_STAR_dummy,
    output:
        "workflow/data/{user}/{project}/alignments/{library}/{library}_Aligned.sortedByCoord.out.bam"
    conda:
        "../envs/master.yaml"
    threads:
        # STAR sometimes fails because of too many opened files (due to high thread count). Lower thread number here if necessary.
        16
    shell:
        """
        STAR --runThreadN {threads} \
             --genomeLoad LoadAndKeep \
             --genomeDir {input.genomeDir} \
             --readFilesIn {input.extracted_fq} \
             --readFilesCommand zcat \
             --limitBAMsortRAM 20000000000 \
             --outFilterMultimapNmax 1 \
             --outFilterType BySJout \
             --outSAMstrandField intronMotif \
             --outFilterIntronMotifs RemoveNoncanonical \
             --outFilterMismatchNmax 6 \
             --outSAMtype BAM SortedByCoordinate \
             --outSAMunmapped Within \
             --outFileNamePrefix workflow/data/{wildcards.user}/{wildcards.project}/alignments/{wildcards.library}/{wildcards.library}_
        """

# Step 4-2: Unload STAR genome index
rule STAR_unload:
    input:
        bams=parse_dynamic_output('STAR'),
        genomeDir=config["genome_index"]
    output:
        temp(touch("tmp/STARunload.done")),
    conda:
        "../envs/master.yaml"
    shell:
        """
        STAR --genomeLoad Remove \
             --genomeDir {input.genomeDir} \
             --outSAMmode None
        rm Log.final.out Log.out Log.progress.out SJ.out.tab
        """

## TODO: STAR mapping summarize statistics

# Step 5-1: Assign reads to genes (featureCount)
rule featurecount:
    input:
        gtf=config["gtf_annotation"],
        bam="workflow/data/{user}/{project}/alignments/{library}/{library}_Aligned.sortedByCoord.out.bam",
        dummy=parse_fc_dummy,
    output:
        assigned="workflow/data/{user}/{project}/alignments/{library}/{library}_gene_assigned",
        bam_counted=temp("workflow/data/{user}/{project}/alignments/{library}/{library}_Aligned.sortedByCoord.out.bam.featureCounts.bam")
    conda:
        "../envs/master.yaml"
    threads:
        config["threads"]
    shell:
        """
        featureCounts -s 1 \
                      -a {input.gtf} \
                      -o {output.assigned} \
                      -R BAM {input.bam} \
                      -T {threads}
        """

# Step 5-2: Assign reads to genes (sort bam files)
rule sambamba_sort:
    input:
        "workflow/data/{user}/{project}/alignments/{library}/{library}_Aligned.sortedByCoord.out.bam.featureCounts.bam"
    output:
        temp("workflow/data/{user}/{project}/alignments/{library}/{library}_assigned_sorted.bam")
    conda:
        "../envs/master.yaml"
    threads:
        config["threads"]
    shell:
        """
        sambamba sort -t {threads} -m 64G \
                      -o {output} \
                      {input}
        """

# Step 6: Count UMIs per gene per cell
rule umi_tools_count:
    input:
        "workflow/data/{user}/{project}/alignments/{library}/{library}_assigned_sorted.bam"
    output:
        temp("workflow/data/{user}/{project}/alignments/{library}/{library}_counts_raw.tsv.gz")
    conda:
        "../envs/master.yaml"
    threads:
        1
    shell:
        """
        umi_tools count --per-gene \
                        --gene-tag=XT \
                        --per-cell \
                        --stdin={input} \
                        --stdout={output}
        """

# Step 7: Append suffix to cells
rule append_sfx:
    input:
        "workflow/data/{user}/{project}/alignments/{library}/{library}_counts_raw.tsv.gz"
    output:
        "workflow/data/{user}/{project}/alignments/{library}/{library}_counts.tsv.gz"
    threads:
        1
    script:
        "../scripts/append_suffix.py"

# Step 8: Aggregate counts
rule aggr_counts:
    input:
        get_files("append_sfx")
    output:
        "workflow/data/{user}/{project}/outs/{project}_counts_all.tsv.gz"
    threads:
        1
    shell:
        """
        zcat {input} | sed '2, ${{/gene/d;}}' | gzip > {output}
        """
