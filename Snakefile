#
configfile: "config.yaml"

rule all:
	input:
                expand("result/fastqc/{sample}_R1_001_fastqc.html", sample=config["samples"]["set2"]),
                expand("result/fastqc/{sample}_R2_001_fastqc.html", sample=config["samples"]["set2"]),
		expand("result/report/stats/{sample}_align.txt", sample=config["samples"]["set2"]),
		expand("result/report/stats/{sample}_markdup.txt", sample=config["samples"]["set2"]),
		expand("result/report/depth/{sample}.mosdepth.summary.txt", sample=config["samples"]["set2"]),
		expand("result/variant/{sample}_autosomes.g.vcf.gz",sample=config["samples"]["set2"])

rule fastqc:
        input:
                "files/{sample}_R1_001.fastq",
                "files/{sample}_R2_001.fastq"
        output:
                "result/fastqc/{sample}_R1_001_fastqc.html",
                "result/fastqc/{sample}_R1_001_fastqc.zip",
                "result/fastqc/{sample}_R2_001_fastqc.html",
                "result/fastqc/{sample}_R2_001_fastqc.zip"
        shell:
                "fastqc -o result/fastqc {input}"

rule bwa_map:
	input:
        	idx=multiext("/userdata/safal/data/ref_genome/masked/GRCh38.p13.genome_masked.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
		read1="/userdata/safal/tutorial/pipeline/files/{sample}_R1_001.fastq",
		read2="/userdata/safal/tutorial/pipeline/files/{sample}_R2_001.fastq"
	output:
		temp("result/alignment/{sample}_mapped.bam")
	params:
		rg=r"@RG\tID:FLOWCELL1.LANE1\tPL:ILLUMINA\tLB:Batch1\tSM:{sample}",
		idx=lambda w, input: os.path.splitext(input.idx[0])[0]
	threads: 28
	log:
		"logs/bwa_mem/{sample}.log"
	shell:
		"bwa mem {params.idx} -R '{params.rg}' -t {threads} {input.read1} {input.read2} | samtools view -Sb - > {output} 2>{log}"


rule samtools_sort:
	input:
		"result/alignment/{sample}_mapped.bam"
	output:
		"result/alignment/{sample}.bam"
	threads: 4
	log:
		"logs/samtools_sort/{sample}.log"
	shell:
		"samtools sort -@ {threads} -T alignment/{wildcards.sample}  "
		"-O bam {input} > {output} 2>{log}"

rule samtools_index:
	input:
		"result/alignment/{sample}.bam"
	output:
		"result/alignment/{sample}.bam.bai"
	threads: 4
	log:
		"logs/samtools_index/{sample}.log"
	shell:
		"samtools index -@ {threads} {input} 2>{log}"

rule mark_dup:
	input:
		bam="result/alignment/{sample}.bam",
		bai="result/alignment/{sample}.bam.bai"
	output:
		out="result/mark_duplicates/{sample}.bam",
		met="result/mark_duplicates/{sample}_metrics.txt"
	log:
		"logs/mark_dup/{sample}.log"
	shell:
		"java -Xmx50g -Djava.io.tmpdir=tmp/{wildcards.sample} -jar /software/picard-2.20.3/picard.jar MarkDuplicates CREATE_INDEX=true I={input.bam} O={output.out} M={output.met} 2>{log}"

rule stats:
        input:
                aln="result/alignment/{sample}.bam",
                dup="result/mark_duplicates/{sample}.bam"
        output:
                aln_out="result/report/stats/{sample}_align.txt",
                dup_out="result/report/stats/{sample}_markdup.txt"
        threads: 4
        shell:
                """
                samtools flagstat -@ {threads} {input.aln} > {output.aln_out}
                samtools flagstat -@ {threads} {input.dup} > {output.dup_out}
                """

rule depth:
        input:
                "result/alignment/{sample}.bam"
        output:
                "result/report/depth/{sample}.mosdepth.summary.txt"
        threads: 4
        shell:
                "mosdepth -t {threads} result/report/depth/{wildcards.sample} {input}"


rule bqsr:
	input:
		ip="result/mark_duplicates/{sample}.bam",
		snp_ref=multiext("/userdata/safal/data/ref_genome/ref_dbsnp.vcf", ".idx"),
		ref=multiext("/userdata/safal/data/ref_genome/masked/GRCh38.p13.genome_masked.fa", ".amb", ".ann", ".bwt", ".pac", ".sa")
	output:
		"result/bqsr/{sample}_recal.table"
	params:
		ref=lambda w, input: os.path.splitext(input.ref[0])[0],
		snp_ref=lambda x, input: os.path.splitext(input.snp_ref[0])[0]
	log:
		"logs/bqsr/{sample}.log"
	shell:
		"/software/gatk-4.1.2.0/gatk --java-options '-Xmx50G' BaseRecalibrator -R {params.ref} --known-sites {params.snp_ref} --tmp-dir /tmp/ -I {input.ip} -O {output} 2>{log}"

rule apply_bqsr:
	input:
		ref=multiext("/userdata/safal/data/ref_genome/masked/GRCh38.p13.genome_masked.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
		file="result/mark_duplicates/{sample}.bam",
		table="result/bqsr/{sample}_recal.table"
	output:
		"result/bqsr/{sample}_bqsr_output.bam"
	params:
		ref=lambda w, input: os.path.splitext(input.ref[0])[0]
	log:
		"logs/apply_bqsr/{sample}.log"
	shell:
		"/software/gatk-4.1.2.0/gatk --java-options '-Xmx50G' ApplyBQSR -R {params.ref} -I {input.file} -OBI true --bqsr-recal-file {input.table} -O {output} 2>{log}"

rule haplotype_caller:
	input:
		ref=multiext("/userdata/safal/data/ref_genome/masked/GRCh38.p13.genome_masked.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
		file2 = "result/bqsr/{sample}_bqsr_output.bam",
		int = "/userdata/safal/sample_SR/interval.bed"
	output:
		"result/variant/{sample}_autosomes.g.vcf.gz"
	params:
		ref=lambda w, input: os.path.splitext(input.ref[0])[0]
	log:
		"logs/haplotype_caller/{sample}.log"
	shell:
		"/software/gatk-4.1.2.0/gatk HaplotypeCaller -R {params.ref} -L {input.int} -I {input.file2} -O {output} -ERC GVCF 2>{log}"

