rule index_genome:
	input:
		os.path.join(config['results'], f'3.Purging/{prefix}.purged.fa')
	output:
		fai=os.path.join(config['results'], f'3.Purging/{prefix}.purged.fa') + ".fai",
		genome=os.path.join(config['results'], f'3.Purging/{prefix}.purged.fa') + ".genome",
		bwt=os.path.join(config['results'], f'3.Purging/{prefix}.purged.fa') + ".bwt",
		pac=os.path.join(config['results'], f'3.Purging/{prefix}.purged.fa') + ".pac",
		ann=os.path.join(config['results'], f'3.Purging/{prefix}.purged.fa') + ".ann",
		amb=os.path.join(config['results'], f'3.Purging/{prefix}.purged.fa') + ".amb",
		sa=os.path.join(config['results'], f'3.Purging/{prefix}.purged.fa') + ".sa"
	conda:
		"yahs"
	threads:
		1
	shell:
		"samtools faidx {input} && cut -f1,2 {input}.fai > {input}.genome && bwa index {input}"

rule map_hic:
	input:
		bwt=os.path.join(config['results'], f'3.Purging/{prefix}.purged.fa') + ".bwt",
		assembly=os.path.join(config['results'], f'3.Purging/{prefix}.purged.fa'),
		reads_r1=expand(hic_R1),
		reads_r2=expand(hic_R2)
	output:
		os.path.join(config['results'], f'4.Scaffolding/{prefix}_bwa.sam')
	conda:
		"yahs"
	threads:
		4
	shell:
		"echo {input.reads_r1} && echo {input.reads_r2} && bwa mem -5SP -T0 -t {threads} -o {output} {input.assembly} <(cat {input.reads_r1}) <(cat {input.reads_r2})" 

rule parse_pairsam: 
	input:
		input_sam=os.path.join(config['results'], f'4.Scaffolding/{prefix}_bwa.sam'),
		input_genome=os.path.join(config['results'], f'3.Purging/{prefix}.purged.fa') + ".genome"
	output:
		os.path.join(config['results'], f'4.Scaffolding/{prefix}_bwa.parsed.pairsam')
	conda:
		"yahs"
	threads:
		4
	shell:
		"pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in {threads} --nproc-out {threads} --chroms-path {input.input_genome} {input.input_sam} > {output}"

rule sort_pairsam:
	input:
		os.path.join(config['results'], f'4.Scaffolding/{prefix}_bwa.parsed.pairsam')
	output:
		os.path.join(config['results'], f'4.Scaffolding/{prefix}_bwa.sorted.pairsam')
	params: 
		tmpdir=os.path.join(config['results'], f'4.Scaffolding')
	conda:
		"yahs"
	threads:
		4	
	shell:
		"pairtools sort --nproc {threads} --tmpdir={params.tmpdir} {input} > {output}"

rule dedup_pairsam:
	input:
		os.path.join(config['results'], f'4.Scaffolding/{prefix}_bwa.sorted.pairsam')
	output:
		output_pairsam=os.path.join(config['results'], f'4.Scaffolding/{prefix}_bwa.dedup.pairsam'),
		output_stats=os.path.join(config['results'], f'4.Scaffolding/{prefix}.dedup.stats')
	conda:
		"yahs"
	threads:
		4
	shell:
		"pairtools dedup --nproc-in {threads} --nproc-out {threads} --mark-dups --output-stats {output.output_stats} --output {output.output_pairsam} {input}"

rule split_pairsam:
	input:
		output_pairsam=os.path.join(config['results'], f'4.Scaffolding/{prefix}_bwa.dedup.pairsam')
	output:
		output_pairs=os.path.join(config['results'], f'4.Scaffolding/{prefix}_bwa.dedup.pairs'),
		output_bam=os.path.join(config['results'], f'4.Scaffolding/{prefix}_bwa.dedup.bam')
	conda:
		"yahs"
	threads:
		4
	shell:
		"pairtools split --nproc-in {threads} --nproc-out {threads} --output-pairs {output.output_pairs} --output-sam {output.output_bam} {input}"

rule sort_bam:
	input:
		os.path.join(config['results'], f'4.Scaffolding/{prefix}_bwa.dedup.bam')
	output:
		os.path.join(config['results'], f'4.Scaffolding/{prefix}_bwa.dedup.sortname.bam')
	params: 
		tmpdir=os.path.join(config['results'], f'4.Scaffolding')
	conda:
		"yahs"
	threads:
		4
	shell:
		"samtools sort -@ {threads} -n -T {params.tmpdir} -o {output} {input}"

rule run_yahs:
	input:
		bam=os.path.join(config['results'], f'4.Scaffolding/{prefix}_bwa.dedup.sortname.bam'),
		assembly=os.path.join(config['results'], f'3.Purging/{prefix}.purged.fa')
	output:
		os.path.join(config['results'], f'4.Scaffolding/{prefix}_scaffolds_final.fa')
	conda:
		"yahs"
	threads:
		1
	shell:
		"yahs {input.assembly} {input.bam} -v 1 -o {prefix}"

