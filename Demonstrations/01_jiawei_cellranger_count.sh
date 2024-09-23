cellranger count \
	--id=cellranger_output_SRR9264343 \
	--transcriptome=references/cellranger_index_v1 \
	--fastqs=fastq/ \
	--sample=SRR9264343 \
	--localcores=8 \
	--localmem=24 \
	--create-bam=true