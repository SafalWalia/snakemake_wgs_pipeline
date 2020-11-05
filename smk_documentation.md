# Documentation
### 1. FastQC:

fastqc -o result/fastqc {input}

* fastqc used from conda environment
* -o result/fastqc : output folder name (both html and zip files would be saved)
* {input}: read file path

Eg: 
```fastqc -o result/fastqc files/ILSHG000001_S1_R1_001.fastq.gz```

### 2. Alignment:
bwa mem {params.idx} -R '{params.rg}' -t {threads} {input.read1} {input.read2} | samtools view -Sb - > {output} 2>{log}

* bwa mem and samtools used from conda environment
* {params.idx}: reference file, with its indexed files
* {params.rg}: read group name inclusion
* {threads}: no of threads used, here 28
* {input.read1} {input.read2}: read file names (R1 & R2)
* {output}: output folder path - the bam file obtained is temprorary (removed once next rule is executed)
* 2>{log}: tool specific log file path, different for each sample  (2> is syntax used)

Eg: 
```bwa mem GRCh38.p13.genome_masked.fa -R "@RG\tID:FLOWCELL1.LANE1\tPL:ILLUMINA\tLB:Batch1\tSM:ILSHG000001_S1" -t 28 ILSHG000001_S1_R1_001.fastq.gz ILSHG000001_S1_R2_001.fastq.gz | samtools view -Sb - > result/alignment/ILSHG000001_S1_mapped.bam 2>logs/bwa_mem/ILSHG000001_S1.log```

### 3. Sort:
samtools sort -@ {threads} -T tmp/alignment/{wildcards.sample} -O bam {input} > {output} 2>{log}

* samtools used from conda environment
* {threads}: no of threads used, here 4
* -T : temp files folder, specific for each sample by its name (hence used - {wildcards.sample})
* {input}, {output}: input and output file names respectively (alignment files) with specified folders

Eg: 
```samtools sort -@ 4 -T tmp/alignment/ILSHG000001_S1 -O bam result/alignment/ILSHG000001_S1_mapped.bam > result/alignment/ILSHG000001_S1.bam 2> logs/samtools_sort/{sample}.log```

### 4. Index:
samtools index -@ {threads} {input} 2>{log}

* samtools used from conda environment
* {threads}: no of threads used, here 4
* {input}: input file name, same as sorted file name with .bai extension

Eg: 
```samtools index -@ 4 tmp/alignment/ILSHG000001_S1.bam 2>logs/samtools_index/ILSHG000001_S1.log```


### 5. Mark duplicates:
java -Xmx50g -Djava.io.tmpdir=tmp/mark_dup/{wildcards.sample} -jar /software/picard-2.20.3/picard.jar MarkDuplicates CREATE_INDEX=true I={input.bam} O={output.out} M={output.met} 2>{log}

* -Xmx50g: max memory allocation
* -Djava.io.tmpdir: to specify temprorary files (name specific)
* picard.jar: used from the software folder
* CREATE_INDEX=true: to generate indexed files
* {input.bam}, {output.out}: input and sorted output bam files path respectively
* {output.met}: metrics file indicating the numbers of duplicates for paired-end reads

Eg: 
```java -Xmx50g -Djava.io.tmpdir=tmp/mark_dup/ILSHG000001_S1 -jar /software/picard-2.20.3/picard.jar MarkDuplicates CREATE_INDEX=true I=result/alignment/ILSHG000001_S1.bam O=result/mark_duplicates/ILSHG000001_S1.bam M=result/mark_duplicates/ILSHG000001_S1_metrics.txt" 2>logs/mark_dup/ILSHG000001_S1.log```


### 6. Alignment Stats:
samtools flagstat -@ {threads} {input.aln} > {output.aln_out}
samtools flagstat -@ {threads} {input.dup} > {output.dup_out}

* samtools used from conda environment
* {threads}: no of threads used, here 4
* {input.aln}, {input.dup}: input file paths before (before marl duplicates) and after duplication respectively
* {output.aln_out}, {output.dup_out}: stats (.txt) file paths for each type of output

Eg:
```samtools flagstat -@ 4 result/alignment/ILSHG000001_S1.bam > result/report/stats/ILSHG000001_S1_align.txt```
```samtools flagstat -@ 4 result/mark_duplicates/ILSHG000001_S1.bam > result/report/stats/ILSHG000001_S1_markdup.txt```


### 7. Mosdepth:
mosdepth -t {threads} result/report/depth/{wildcards.sample} {input}

* mosdepth used from local folder 
* {threads}: no of threads used, here 4
* result/report/depth/{wildcards.sample}: output folder, name specific 
* {input}: path of input files obtained after alignment(sorted and indexed)

Eg:
```mosdepth -t 4 result/report/depth/ILSHG000001_S1 result/alignment/ILSHG000001_S1.bam```



### 8. BQSR:
/software/gatk-4.1.2.0/gatk --java-options '-Xmx50G' BaseRecalibrator -R {params.ref} --known-sites {params.snp_ref} --tmp-dir /tmp/ -I {input.ip} -O {output} 2>{log}

* gatk used from software folder
* {params.ref}: reference genome path (with all extensions available after indexing)
* {params.snp_ref}: indexed reference variant file from dbSNP 
* --tmp-dir: path of temprorary files directory
* {input.ip}: input file path (bam file obtained after marked duplicates)
* {output}: path for obtained recalibration table 

Eg:
```/software/gatk-4.1.2.0/gatk --java-options '-Xmx50G' BaseRecalibrator -R GRCh38.p13.genome_masked.fa  --known-sites ref_dbsnp.vcf ref_dbsnp.vcf.idx --tmp-dir tmp/bqsr -I result/mark_duplicates/ILSHG000001_S1.bam -O result/bqsr/ILSHG000001_S1_recal.table 2>logs/bqsr/ILSHG000001_S1.log```



### 9. Apply BQSR:
/software/gatk-4.1.2.0/gatk --java-options '-Xmx50G' ApplyBQSR -R {params.ref} -I {input.file} -OBI true --bqsr-recal-file {input.table} -O {output} 2>{log}

* gatk used from software folder
* {params.ref}: reference genome path (with all extensions available after indexing)
* {input.file}: input file
* -OBI true: to enable indexing
* --bqsr-recal-file {input.table}: path of recalibration file
* {output}: sorted BAM file 
* {log}: log file

Eg:
```/software/gatk-4.1.2.0/gatk --java-options '-Xmx50G' ApplyBQSR -R GRCh38.p13.genome_masked.fa -I result/mark_duplicates/ILSHG000001_S1.bam -OBI true --bqsr-recal-file result/bqsr/ILSHG000001_S1_recal.table -O result/bqsr/ILSHG000001_S1_bqsr_output.bam 2>logs/apply_bqsr/ILSHG000001_S1.log```

### 10. Haplotype Caller:
/software/gatk-4.1.2.0/gatk HaplotypeCaller -R {params.ref} -L {input.int} -I {input.file2} -O {output} -ERC GVCF 2>{log}

* gatk used from software folder
* {params.ref}: reference genome path (with all extensions available after indexing)
* {input.int}: Interval file in BED format, with co-ordinates of autosomes only
* {input.file2}: Input BAM file after BQSR
* {output}: compressed GVCF file
* {log}: log file

Eg:
```/software/gatk-4.1.2.0/gatk HaplotypeCaller -R GRCh38.p13.genome_masked.fa -L interval.bed -I result/bqsr/ILSHG000001_S1_bqsr_output.bam -O result/variant/ILSHG000001_S1_autosomes.g.vcf.gz -ERC GVCF 2>logs/haplotype_caller/ILSHG000001_S1.log```
