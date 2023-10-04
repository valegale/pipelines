def contig_output(contig_assembler):
	"""
	Three possible assembler: hifiasm, flye, hicanu. 
	"""
	if contig_assembler == "hifiasm":
		asm_p_fa = os.path.join(config['results'], '2.Contigging/hifiasm/purgeL' + config['purgel'] + '/hifiasmL' + config['purgel'] + '.asm.p_ctg.fasta')
		return asm_p_fa
	elif contig_assembler == "flye":
		asm_flye = os.path.join(config['results'], '2.Contigging/flye/assembly.fasta')
		return asm_flye
	elif contig_assembler == "hicanu":
		asm_hicanu = os.path.join(config['results'], '2.Contigging/hicanu/asm.contigs.fasta')
		return asm_hicanu

rule split_ref:
	input:
		contig_output(contig_assembler)
	output:
		os.path.join(config['results'], f'3.Purging/{prefix}_split.fasta')
	conda:
		os.path.join(config['snakemake_dir_path'], 'envs/purging.yaml')
	shell:
		"split_fa {input} > {output}"

rule self_map:
	input:
		os.path.join(config['results'], f'3.Purging/{prefix}_split.fasta')
	output:
		os.path.join(config['results'], f'3.Purging/{prefix}_split.paf')
	conda:
		os.path.join(config['snakemake_dir_path'], 'envs/purging.yaml')
	threads: 16
	shell:
		"minimap2 -I200G -t {threads} -xasm5 -DP {input} {input} > {output}"

rule map_reads:
	input:
		ref=contig_output(contig_assembler),
		reads=lambda wildcards: samples[samples["sample_id"] == wildcards.sample_id]["file_path"].values[0]
	output:
		os.path.join(config['results'], '3.Purging/{sample_id}.paf')
	conda:
		os.path.join(config['snakemake_dir_path'], 'envs/purging.yaml')
	threads: 16
	shell:
		"minimap2 -I 200G -x map-pb -t {threads} {input.ref} {input.reads} > {output}"

rule pbstats:
	input:
		expand(os.path.join(config['results'], '3.Purging/{sample_id}.paf'), sample_id=samples["sample_id"])
	output:
		os.path.join(config['results'], f'3.Purging/{prefix}_coverage/PB.stat'),
		os.path.join(config['results'], f'3.Purging/{prefix}_coverage/PB.cov.wig'),
		os.path.join(config['results'], f'3.Purging/{prefix}_coverage/PB.base.cov')
	params:
		outputDir=f'{prefix}_coverage'
	conda:
		os.path.join(config['snakemake_dir_path'], 'envs/purging.yaml')
	shell:
		"pbcstat -O {params.outputDir} {input}"

rule calcuts:
	input:
		os.path.join(config['results'], f'3.Purging/{prefix}_coverage/PB.stat')
	output:
		os.path.join(config['results'], f'3.Purging/{prefix}_cutoffs'),
	conda:
		os.path.join(config['snakemake_dir_path'], 'envs/purging.yaml')
	shell:
		"calcuts {input} > {output}"

rule purge_dups:
	input:
		cutoffs=os.path.join(config['results'], f'3.Purging/{prefix}_cutoffs'),
		coverage=os.path.join(config['results'], f'3.Purging/{prefix}_coverage/PB.base.cov'),
		split_paf=os.path.join(config['results'], f'3.Purging/{prefix}_split.paf')
	output:
		os.path.join(config['results'], f'3.Purging/{prefix}_dups.bed')
	conda:
		os.path.join(config['snakemake_dir_path'], 'envs/purging.yaml')
	shell:
		"purge_dups -2 -c {input.coverage} -T {input.cutoffs} {input.split_paf} > {output}"

rule get_seqs:
	input:
		dups=os.path.join(config['results'], f'3.Purging/{prefix}_dups.bed'),
		genome=contig_output(contig_assembler)
	output:
		os.path.join(config['results'], f'3.Purging/{prefix}.purged.fa'),
		os.path.join(config['results'], f'3.Purging/{prefix}.hap.fa')
	conda:
		os.path.join(config['snakemake_dir_path'], 'envs/purging.yaml')
	shell:
		"get_seqs -e -p {prefix} {input.dups} {input.genome}"

