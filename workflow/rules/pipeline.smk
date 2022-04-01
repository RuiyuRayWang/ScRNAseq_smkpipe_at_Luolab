# Step 0: Aggregate fastqs
# rule aggr_fqs:
#     input:
#         get_table
#     output:
#         "workflow/data/{user}/{project}/fastqs/{sample}/{sample}_R1.fq.gz",
#         "workflow/data/{user}/{project}/fastqs/{sample}/{sample}_R2.fq.gz",
#     threads:
#         1
#     script:
#         "../scripts/aggr_fqs.py"

rule aggr_r1_fqs:
    input:
        unpack(get_r1_aggr_input)
    output:
        "workflow/data/{user}/{project}/fastqs/{sample}/{sample}_R1.fq.gz",
    threads:
        1
    shell:
        "cat {input} > {output}"

rule aggr_r2_fqs:
    input:
        unpack(get_r2_aggr_input)
    output:
        "workflow/data/{user}/{project}/fastqs/{sample}/{sample}_R2.fq.gz",
    threads:
        1
    shell:
        "cat {input} > {output}"

# Step 1: Identify cell barcode whitelist (identify correct BC)
rule umi_tools_whitelist:
    input:
        unpack(get_fastqs),
        # tmp="tmp/aggr_fqs.done",
    output:
        "workflow/data/{user}/{project}/alignments/{sample}/{sample}_whitelist.txt",
        report("workflow/data/{user}/{project}/alignments/{sample}/{sample}_cell_barcode_counts.png", caption="../report/umi_tools_whitelist.rst", category="umi_tools")
    log:
        "workflow/data/{user}/{project}/logs/{sample}/{sample}_whitelist.log"
    conda:
        "../envs/master.yaml"
    threads:
        1
    shell:
        """
        umi_tools whitelist --bc-pattern=CCCCCCCCNNNNNNNN \
                            --log={log} \
                            --set-cell-number=100 \
                            --plot-prefix=workflow/data/{wildcards.user}/{wildcards.project}/alignments/{wildcards.sample}/{wildcards.sample} \
                            --stdin={input.read1} \
                            --stdout={output}
        """

# Step 2: Wash whitelist
rule wash_whitelist:
    input:
        whitelist="workflow/data/{user}/{project}/alignments/{sample}/{sample}_whitelist.txt",
        ground_truth=config['ground_truth']
    output:
        "workflow/data/{user}/{project}/alignments/{sample}/{sample}_whitelist_washed.txt"
    threads:
        1
    script:
        "../scripts/wash_whitelist.py"

# Step 3: Extract barcodes and UMIs and add to read names
rule umi_tools_extract:
    input:
        unpack(get_fastqs),
        whitelist_washed="workflow/data/{user}/{project}/alignments/{sample}/{sample}_whitelist_washed.txt"
    output:
        "workflow/data/{user}/{project}/alignments/{sample}/{sample}_extracted.fq.gz",
    log:
        "workflow/data/{user}/{project}/logs/{sample}/{sample}_extract.log"
    conda:
        "../envs/master.yaml"
    threads:
        1
    shell:
        """
        umi_tools extract --bc-pattern=CCCCCCCCNNNNNNNN \
                          --log={log} \
                          --stdin={input.read1} \
                          --read2-in={input.read2} \
                          --stdout={output} \
                          --read2-stdout \
                          --filter-cell-barcode \
                          --error-correct-cell \
                          --whitelist={input.whitelist_washed}
        """

# Step 4-0: Generate genome index
rule STAR_gen:
    input:
        genome_fa=config["genome_fa"],
        gtf=config["gtf_annotation"]
    output:
        directory(config["genome_index"])
    params:
        sjdbOverhang=config["sjdbOverhang"]
    cache: True
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

# Step 4-1: Load genome index to memory
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

# Step 4-2: Map reads
rule STAR:
    input:
        extracted_fq="workflow/data/{user}/{project}/alignments/{sample}/{sample}_extracted.fq.gz",
        genomeDir=config["genome_index"],
        # dummy="tmp/STARload.done"
        dummy=parse_STAR_dummy,
    output:
        bam="workflow/data/{user}/{project}/alignments/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
        lf=report("workflow/data/{user}/{project}/alignments/{sample}/{sample}_Log.final.out", caption="../report/STAR.rst", category="STAR")
    conda:
        "../envs/master.yaml"
    threads:
        # STAR sometimes fails as too many files are opened (due to high thread number). Decrease thread number here if necessary.
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
             --outFileNamePrefix workflow/data/{wildcards.user}/{wildcards.project}/alignments/{wildcards.sample}/{wildcards.sample}_
        """

# Step 4-3: Unload STAR genome
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
        """

# Step 5: Assign reads to genes (featureCount)
rule featureCounts:
    input:
        gtf=config["gtf_annotation"],
        bam="workflow/data/{user}/{project}/alignments/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
        # dummy="tmp/STARunload.done"
        dummy=parse_fc_dummy,
    output:
        assigned=temp("workflow/data/{user}/{project}/alignments/{sample}/{sample}_gene_assigned"),
        summary=report("workflow/data/{user}/{project}/alignments/{sample}/{sample}_gene_assigned.summary", caption="../report/featureCounts.rst", category="featureCounts"),
        bam_counted=temp("workflow/data/{user}/{project}/alignments/{sample}/{sample}_Aligned.sortedByCoord.out.bam.featureCounts.bam")
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

# Step 6: Sort and index BAM
rule sambamba_sort:
    input:
        "workflow/data/{user}/{project}/alignments/{sample}/{sample}_Aligned.sortedByCoord.out.bam.featureCounts.bam"
    output:
        bam=temp("workflow/data/{user}/{project}/alignments/{sample}/{sample}_assigned_sorted.bam"),
        index=temp("workflow/data/{user}/{project}/alignments/{sample}/{sample}_assigned_sorted.bam.bai")
    conda:
        "../envs/master.yaml"
    threads:
        config["threads"]
    shell:
        """
        sambamba sort -t {threads} -m 64G \
                      -o {output.bam} \
                      {input}
        """

# (Optional) Step 7-1: Parse CB UB tag
rule pysam:
    input:
        "workflow/data/{user}/{project}/alignments/{sample}/{sample}_assigned_sorted.bam",
        "workflow/data/{user}/{project}/alignments/{sample}/{sample}_assigned_sorted.bam.bai"
    output:
        temp("workflow/data/{user}/{project}/alignments/{sample}/{sample}_tagged.bam")
    conda:
        "../envs/velocyto.yaml"
    threads:
        1
    script:
        "../scripts/edit_bam_tag.py"

# Step 7-2: Create pre-sorted BAM for velocyto
rule pre_sort:
    input:
        "workflow/data/{user}/{project}/alignments/{sample}/{sample}_tagged.bam"
    output:
        temp("workflow/data/{user}/{project}/alignments/{sample}/cellsorted_{sample}_tagged.bam")
    conda:
        "../envs/master.yaml"
    threads:
        config["threads"]
    shell:
        """
        samtools sort -@ {threads} \
                      -t CB \
                      -O BAM \
                      -o {output} \
                      {input}
        """

# Step 7-3: Generate velocyto loom file
rule velocyto:
    input:
        bam="workflow/data/{user}/{project}/alignments/{sample}/{sample}_tagged.bam",
        sorted_bam="workflow/data/{user}/{project}/alignments/{sample}/cellsorted_{sample}_tagged.bam",
        gtf=config["gtf_annotation"]
    output:
        temp("workflow/data/{user}/{project}/alignments/{sample}/{sample}_velocyto.loom")
    conda:
        "../envs/velocyto.yaml"
    threads:
        16
    shell:
        """
        velocyto run -o workflow/data/{wildcards.user}/{wildcards.project}/alignments/{wildcards.sample} \
                     -v \
                     -e {wildcards.sample}_velocyto \
                     {input.bam} \
                     {input.gtf}
        """

# Step 7-4: Aggregate velocyto loom file
rule aggr_velo_loom:
    input:
        get_files("velocyto")
    output:
        "workflow/data/{user}/{project}/outs/{project}_velocyto_all.loom"
    conda:
        "../envs/velocyto.yaml"
    threads:
        1
    script:
        "../scripts/aggr_velo_loom.py"

# Step 8: Count UMIs per gene per cell
rule umi_tools_count:
    input:
        bam="workflow/data/{user}/{project}/alignments/{sample}/{sample}_assigned_sorted.bam",
        index="workflow/data/{user}/{project}/alignments/{sample}/{sample}_assigned_sorted.bam.bai"
    output:
        temp("workflow/data/{user}/{project}/alignments/{sample}/{sample}_counts_raw.tsv.gz")
    log:
        "workflow/data/{user}/{project}/logs/{sample}/{sample}_count.log"
    conda:
        "../envs/master.yaml"
    threads:
        1
    shell:
        """
        umi_tools count --per-gene \
                        --per-cell \
                        --gene-tag=XT \
                        --log={log} \
                        --stdin={input.bam} \
                        --stdout={output} \
        """

# Step 9: Append suffix to cells
rule append_sfx:
    input:
        "workflow/data/{user}/{project}/alignments/{sample}/{sample}_counts_raw.tsv.gz"
    output:
        "workflow/data/{user}/{project}/alignments/{sample}/{sample}_counts.tsv.gz"
    threads:
        1
    script:
        "../scripts/append_suffix.py"

# Step 10: Aggregate counts
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

# Step 11: QC report
rule qc_report:
    input:
        log_whitelist=get_logfiles("log_whitelist"),
        log_extract=get_logfiles("log_extract"),
        log_count=get_logfiles("log_count")
    output:
        df_stats=report("workflow/data/{user}/{project}/outs/{project}_stats.csv", caption="../report/stats_table.rst", category="Aggregated Stats"),
        reads_pdf=report("workflow/data/{user}/{project}/outs/{project}_reads_stats.pdf", caption="../report/reads_stats.rst", category="Aggregated Stats"),
        reads_html="workflow/data/{user}/{project}/outs/{project}_reads_stats.html",
    conda:
        "../envs/master.yaml"
    threads:
        1
    script:
        "../scripts/parse_qc_stats.py"

onsuccess:
    shell("rm -f Log.final.out Log.out Log.progress.out SJ.out.tab geckodriver.log")

# onerror:
#     ## TODO: Get conda env through parser
#     shell(
#         "rm -f Log.final.out Log.out Log.progress.out SJ.out.tab geckodriver.log" + \
#         "STAR --genomeLoad Remove --genomeDir {config[genome_index]} --outSAMmode None;",
#         # "mail -s ", 
#         # "STAR --version",
#         conda_env="/home/luolab/GITHUB_REPO/ScRNAseq_smkpipe_at_Luolab/.snakemake/conda/e9e317d19678f2458891b4f29ae2546c")
