import csv
import os
from collections import OrderedDict
from pathlib import Path


configfile: os.environ.get("BULK_RNASEQ_CONFIG", "config/config.yaml")


def cfg(name, default=None):
    return config.get(name, default)


def clean_value(value):
    return (value or "").strip()


def clean_path(value):
    return clean_value(value).replace("\\", "/")


SAMPLE_TABLE = clean_path(cfg("sample_table", "config/samples.tsv"))
RUN_NAME = clean_value(cfg("run_name", Path(SAMPLE_TABLE).stem))
THREADS = int(cfg("threads", 8))
STRANDEDNESS = int(cfg("strandedness", 0))

DIRS = cfg("directories", {})
TRIM_DIR = clean_path(DIRS.get("trimmed_fastq", "data/trimmed_fastq"))
QC_DIR = clean_path(DIRS.get("qc", "results/qc"))
BAM_DIR = clean_path(DIRS.get("bam", "results/bam"))
COUNTS_DIR = clean_path(DIRS.get("counts", "results/counts"))
DESEQ2_DIR = clean_path(DIRS.get("deseq2", "results/deseq2"))
GO_DIR = clean_path(DIRS.get("go", "results/go"))
GSEA_DIR = clean_path(DIRS.get("gsea", "results/gsea"))
METADATA_DIR = clean_path(DIRS.get("metadata", "results/metadata"))
LOG_DIR = clean_path(DIRS.get("logs", "logs"))

REF = cfg("reference", {})
GENOME_FA = clean_path(REF.get("genome_fa", "data/reference/GRCh38.primary_assembly.genome.fa"))
GTF_FILE = clean_path(REF.get("gtf", "data/reference/gencode.v44.annotation.gtf"))
HISAT2_INDEX = clean_path(REF.get("hisat2_index_prefix", "data/reference/hisat2_index/GRCh38"))

HISAT2_CFG = cfg("hisat2", {})
HISAT2_THREADS = int(HISAT2_CFG.get("threads", 2))
SORT_THREADS = int(HISAT2_CFG.get("sort_threads", 1))
SORT_MEMORY = clean_value(HISAT2_CFG.get("sort_memory", "768M"))

COMPARISONS = cfg("comparisons", [])
COMPARISON_NAMES = [clean_value(c["name"]) for c in COMPARISONS]
COMPARISON_ARG = ",".join(
    [
        f"{clean_value(c['numerator_group'])}:{clean_value(c['denominator_group'])}:{clean_value(c['name'])}"
        for c in COMPARISONS
    ]
)


def load_samples(sample_table):
    required = ["run_id", "cell_line", "time", "treatment", "fq"]
    records = OrderedDict()

    with open(sample_table, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != required:
            raise ValueError(
                "Invalid sample table header. Expected: "
                + "\t".join(required)
                + "; found: "
                + "\t".join(reader.fieldnames or [])
            )

        for row in reader:
            row = {key: clean_value(row.get(key)) for key in required}
            if not row["run_id"]:
                continue
            row["fq"] = clean_path(row["fq"])

            previous = records.get(row["run_id"])
            if previous is not None:
                same_record = all(previous[key] == row[key] for key in required)
                if same_record:
                    continue
                raise ValueError(
                    f"run_id {row['run_id']} appears more than once with different metadata."
                )

            records[row["run_id"]] = row

    if not records:
        raise ValueError(f"No samples found in {sample_table}")

    return records


SAMPLE_RECORDS = load_samples(SAMPLE_TABLE)
SAMPLES = list(SAMPLE_RECORDS.keys())
FQ_BY_SAMPLE = {sample: SAMPLE_RECORDS[sample]["fq"] for sample in SAMPLES}

METADATA_TABLE = f"{METADATA_DIR}/{RUN_NAME}.samples.cleaned.tsv"
INDEX_DONE = f"{HISAT2_INDEX}.index.done"


rule all:
    input:
        f"{QC_DIR}/{RUN_NAME}/multiqc/multiqc_report.html",
        expand(f"{TRIM_DIR}/{RUN_NAME}/{{sample}}.trimmed.fastq.gz", sample=SAMPLES),
        expand(f"{BAM_DIR}/{RUN_NAME}/{{sample}}.sorted.bam", sample=SAMPLES),
        expand(f"{BAM_DIR}/{RUN_NAME}/{{sample}}.sorted.bam.bai", sample=SAMPLES),
        f"{COUNTS_DIR}/{RUN_NAME}/gene_counts.tsv",
        f"{DESEQ2_DIR}/.deseq2.done",
        f"{GO_DIR}/.go.done",
        f"{GSEA_DIR}/.gsea.done"


rule sample_metadata:
    output:
        METADATA_TABLE
    run:
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        with open(output[0], "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["run_id", "cell_line", "time", "treatment", "fq"],
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(SAMPLE_RECORDS.values())


rule fastqc_raw:
    input:
        fq=lambda wildcards: FQ_BY_SAMPLE[wildcards.sample]
    output:
        html=f"{QC_DIR}/{RUN_NAME}/fastqc/{{sample}}_fastqc.html",
        zip=f"{QC_DIR}/{RUN_NAME}/fastqc/{{sample}}_fastqc.zip"
    log:
        f"{LOG_DIR}/01_fastqc/{RUN_NAME}/{{sample}}.log"
    threads: THREADS
    shell:
        """
        mkdir -p {QC_DIR}/{RUN_NAME}/fastqc $(dirname {log})
        fastqc {input.fq} --threads {threads} --outdir {QC_DIR}/{RUN_NAME}/fastqc > {log} 2>&1
        """


rule fastp:
    input:
        fq=lambda wildcards: FQ_BY_SAMPLE[wildcards.sample]
    output:
        fq=f"{TRIM_DIR}/{RUN_NAME}/{{sample}}.trimmed.fastq.gz",
        json=f"{QC_DIR}/{RUN_NAME}/fastp/{{sample}}.fastp.json",
        html=f"{QC_DIR}/{RUN_NAME}/fastp/{{sample}}.fastp.html"
    log:
        f"{LOG_DIR}/02_fastp/{RUN_NAME}/{{sample}}.log"
    threads: THREADS
    shell:
        """
        mkdir -p {TRIM_DIR}/{RUN_NAME} {QC_DIR}/{RUN_NAME}/fastp $(dirname {log})
        fastp --in1 {input.fq} --out1 {output.fq} --thread {threads} --json {output.json} --html {output.html} > {log} 2>&1
        """


rule multiqc:
    input:
        fastqc=expand(f"{QC_DIR}/{RUN_NAME}/fastqc/{{sample}}_fastqc.html", sample=SAMPLES),
        fastp=expand(f"{QC_DIR}/{RUN_NAME}/fastp/{{sample}}.fastp.html", sample=SAMPLES)
    output:
        f"{QC_DIR}/{RUN_NAME}/multiqc/multiqc_report.html"
    log:
        f"{LOG_DIR}/03_multiqc/{RUN_NAME}/multiqc.log"
    shell:
        """
        mkdir -p {QC_DIR}/{RUN_NAME}/multiqc $(dirname {log})
        multiqc {QC_DIR}/{RUN_NAME} --outdir {QC_DIR}/{RUN_NAME}/multiqc --filename multiqc_report.html --force > {log} 2>&1
        """


rule hisat2_index:
    input:
        genome=GENOME_FA
    output:
        touch(INDEX_DONE)
    log:
        f"{LOG_DIR}/04_hisat2_index/hisat2-build.log"
    threads: THREADS
    shell:
        """
        mkdir -p $(dirname {HISAT2_INDEX}) $(dirname {log})
        hisat2-build -p {threads} {input.genome} {HISAT2_INDEX} > {log} 2>&1
        if [ ! -s {HISAT2_INDEX}.1.ht2 ] && [ ! -s {HISAT2_INDEX}.1.ht2l ]; then
            echo "HISAT2 index was not created for prefix {HISAT2_INDEX}" >> {log}
            exit 1
        fi
        """


rule hisat2_align:
    input:
        fq=f"{TRIM_DIR}/{RUN_NAME}/{{sample}}.trimmed.fastq.gz",
        index_done=INDEX_DONE
    output:
        bam=f"{BAM_DIR}/{RUN_NAME}/{{sample}}.sorted.bam",
        bai=f"{BAM_DIR}/{RUN_NAME}/{{sample}}.sorted.bam.bai"
    params:
        index=HISAT2_INDEX,
        sort_memory=SORT_MEMORY
    log:
        align=f"{LOG_DIR}/05_hisat2_align/{RUN_NAME}/{{sample}}.log",
        summary=f"{LOG_DIR}/05_hisat2_align/{RUN_NAME}/{{sample}}.summary.txt"
    threads: HISAT2_THREADS + SORT_THREADS
    shell:
        """
        set -euo pipefail
        mkdir -p {BAM_DIR}/{RUN_NAME} $(dirname {log.align})
        tmp_bam={BAM_DIR}/{RUN_NAME}/{wildcards.sample}.sorted.tmp.bam
        rm -f "$tmp_bam" "$tmp_bam.bai"
        hisat2 -x {params.index} -U {input.fq} -p {HISAT2_THREADS} --summary-file {log.summary} 2> {log.align} | \
            samtools sort -@ {SORT_THREADS} -m {params.sort_memory} -T {BAM_DIR}/{RUN_NAME}/{wildcards.sample}.sort_tmp -o "$tmp_bam" - 2>> {log.align}
        samtools quickcheck "$tmp_bam"
        mv "$tmp_bam" {output.bam}
        samtools index {output.bam} >> {log.align} 2>&1
        """


rule featurecounts:
    input:
        bams=expand(f"{BAM_DIR}/{RUN_NAME}/{{sample}}.sorted.bam", sample=SAMPLES),
        gtf=GTF_FILE
    output:
        counts=f"{COUNTS_DIR}/{RUN_NAME}/gene_counts.tsv",
        summary=f"{COUNTS_DIR}/{RUN_NAME}/gene_counts.tsv.summary"
    log:
        f"{LOG_DIR}/06_featurecounts/{RUN_NAME}/featurecounts.log"
    threads: THREADS
    shell:
        """
        mkdir -p {COUNTS_DIR}/{RUN_NAME} $(dirname {log})
        featureCounts -T {threads} -a {input.gtf} -t exon -g gene_id -s {STRANDEDNESS} -o {output.counts} {input.bams} > {log} 2>&1
        """


rule deseq2:
    input:
        counts=f"{COUNTS_DIR}/{RUN_NAME}/gene_counts.tsv",
        samples=METADATA_TABLE
    output:
        touch(f"{DESEQ2_DIR}/.deseq2.done")
    params:
        out_dir=DESEQ2_DIR,
        comparisons=COMPARISON_ARG
    log:
        f"{LOG_DIR}/07_deseq2/{RUN_NAME}.log"
    shell:
        """
        mkdir -p {params.out_dir} $(dirname {log})
        Rscript scripts/07_deseq2.R {input.counts} {input.samples} {params.out_dir} "{params.comparisons}" > {log} 2>&1
        """


rule go_enrichment:
    input:
        deseq2_done=f"{DESEQ2_DIR}/.deseq2.done"
    output:
        touch(f"{GO_DIR}/.go.done")
    params:
        deseq2_dir=DESEQ2_DIR,
        out_dir=GO_DIR,
        comparisons=COMPARISON_ARG
    log:
        f"{LOG_DIR}/08_go/{RUN_NAME}.log"
    shell:
        """
        mkdir -p {params.out_dir} $(dirname {log})
        Rscript scripts/08_go.R {params.deseq2_dir} {params.out_dir} "{params.comparisons}" > {log} 2>&1
        """


rule gsea:
    input:
        deseq2_done=f"{DESEQ2_DIR}/.deseq2.done",
        go_done=f"{GO_DIR}/.go.done"
    output:
        touch(f"{GSEA_DIR}/.gsea.done")
    params:
        deseq2_dir=DESEQ2_DIR,
        go_dir=GO_DIR,
        out_dir=GSEA_DIR,
        comparisons=COMPARISON_ARG
    log:
        f"{LOG_DIR}/09_gsea/{RUN_NAME}.log"
    shell:
        """
        mkdir -p {params.out_dir} $(dirname {log})
        Rscript scripts/09_gsea.R {params.deseq2_dir} {params.go_dir} {params.out_dir} "{params.comparisons}" > {log} 2>&1
        """
