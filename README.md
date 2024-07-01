# WGS
Analysis code for genome data

## Data processing

### 1.QC

    fastp -w 8 -i sample1/sample1_r1.fq.gz -I sample1/sample1_r2.fq.gz -o QC/sample1/sample1_R1.clean.fq.gz -O QC/sample1/sample1_R2.clean.fq.gz --html QC/sample1/sample1.html --json QC/sample1/sample1.json

### 2.Mapping

    (sentieon bwa mem -t 42 -K 10000000 -M -R "@RG\tID:sample1\tSM:sample1\tLB:PE100\tPU:PE100\tPL:illumina\tCN:WAW" genome.fa sample1/sample1_R1.clean.fq.gz sample1/sample1_R2.clean.fq.gz) | sentieon util sort -r genome.fa -o ./DNAdb/sample1/sample1.sorted.bam -t 42 --sam2bam --intermediate_compress_level 6 -i -

### 3.Metrics

    sentieon driver -r genome.fa -t 42 -i ./DNAdb/sample1/sample1.sorted.bam --algo MeanQualityByCycle ./DNAdb/sample1/metrics/sample1_mq_metrics.txt --algo QualDistribution ./DNAdb/sample1/metrics/sample1_qd_metrics.txt --algo GCBias  --summary ./DNAdb/sample1/metrics/sample1_gc_summary.txt ./DNAdb/sample1/metrics/sample1_gc_metrics.txt --algo AlignmentStat --adapter_seq '' ./DNAdb/sample1/metrics/sample1_aln_metrics.txt --algo InsertSizeMetricAlgo ./DNAdb/sample1/metrics/sample1_is_metrics.txt
    sentieon plot metrics -o ./DNAdb/sample1/metrics/sample1_metrics-report.pdf gc=./DNAdb/sample1/metrics/sample1_gc_metrics.txt qd=./DNAdb/sample1/metrics/sample1_qd_metrics.txt mq=./DNAdb/sample1/metrics/sample1_mq_metrics.txt isize=./DNAdb/sample1/metrics/sample1_is_metrics.txt

### 4.Deduplication

    sentieon driver  -t 42 -i ./DNAdb/sample1/sample1.sorted.bam  --algo LocusCollector --fun score_info ./DNAdb/sample1/metrics/sample1_score.txt
    sentieon driver  -t 42 -i ./DNAdb/sample1/sample1.sorted.bam  --algo Dedup  --bam_compression 6 --rmdup --score_info ./DNAdb/sample1/metrics/sample1_score.txt --metrics ./DNAdb/sample1/metrics/sample1_dedup_metrics.txt ./DNAdb/sample1/sample1.deduped.bam
    sentieon driver -r genome.fa -t 42 -i ./DNAdb/sample1/sample1.deduped.bam --algo CoverageMetrics --omit_interval_stat --omit_sample_stat --omit_base_output ./DNAdb/sample1/metrics/sample1_cov_metrics.txt

### 5.Realigner

    sentieon driver -r genome.fa -t 42 -i ./DNAdb/sample1/sample1.deduped.bam --algo Realigner --bam_compression 6 ./DNAdb/sample1/sample1.realn.bam

### 6.Base recalibration

    sentieon driver -r genome.fa -t 42 -i ./DNAdb/sample1/sample1.realn.bam --algo QualCal ./DNAdb/sample1/metrics/sample1_recal_data.table
    sentieon driver -r genome.fa -t 42 -i ./DNAdb/sample1/sample1.realn.bam  -q ./DNAdb/sample1/metrics/sample1_recal_data.table --algo QualCal ./DNAdb/sample1/metrics/sample1_recal_data.table.post  --algo ReadWriter --bam_compression 8 ./DNAdb/sample1/sample1.recaled.bam
    sentieon driver -t 42 --algo QualCal --plot --before ./DNAdb/sample1/metrics/sample1_recal_data.table --after ./DNAdb/sample1/metrics/sample1_recal_data.table.post ./DNAdb/sample1/metrics/sample1_recal.csv
    sentieon plot bqsr -o ./DNAdb/sample1/metrics/sample1_recal_plots.pdf ./DNAdb/sample1/metrics/sample1_recal.csv

### 7.HC Variant calling

    sentieon driver -r genome.fa -t 42 -i ./DNAdb/sample1/sample1.recaled.bam --algo DNAscope --emit_mode gvcf --var_type snp,indel  ./DNAdb/sample1/sample1_DNAscope.GVCF.gz

### 8.Concat

    sentieon driver  -r genome.fa -t 128 --algo GVCFtyper -v sample1.gvcf.gz -v sample2.gvcf.gz -v sample3.gvcf.gz  gatk.raw.vcf

### 9.Filter

    vcftools --vcf gatk.raw.vcf  --min-meanDP 5 --max-missing-count 300  --recode --recode-INFO-all --out gatk
    gatk --java-options SelectVariants -R genome.fa -V gatk.recode.vcf --select-type-to-include SNP -O snp.raw.vcf
    gatk --java-options SelectVariants -R genome.fa -V gatk.recode.vcf --select-type-to-include INDEL -O indel.raw.vcf
    gatk --java-options VariantFiltration -R genome.fa -V snp.raw.vcf --filter-expression "QUAL < 30.0 || QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 4.0 || ReadPosRankSum < -8.0"  --filter-name "my_snp_filter" --missing-values-evaluate-as-failing true -O snp.raw.tmp.vcf
    gatk --java-options VariantFiltration -R genome.fa -V indel.raw.vcf --filter-expression "QUAL < 30.0 || QD < 2.0 || FS > 200.0 || SOR > 10.0 || ReadPosRankSum < -20.0 || MQ < 40.0 || MQRankSum < -12.5" --filter-name "my_indel_filter" --missing-values-evaluate-as-failing true -O indel.raw.tmp.vcf
    gatk --java-options SelectVariants -R genome.fa -V snp.raw.tmp.vcf --exclude-filtered -O snp.filter.vcf
    gatk --java-options SelectVariants -R genome.fa -V indel.raw.tmp.vcf --exclude-filtered -O indel.filter.vcf

