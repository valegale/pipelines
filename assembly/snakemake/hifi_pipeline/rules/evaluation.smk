rule merqury:
	""" Execute Merqury following Hifiasm to evaluate the quality and accuracy of the assembled genome."""

	input:
		asm_p = contig_output_primary(contig_assembler, phasing_mode),
		meryl_db = os.path.join(config['results'] , prefix, "genome_profiling", (prefix +".meryl"))
	output:
		out_qv = os.path.join(evaluation_folder(contig_assembler, phasing_mode), 'merqury.qv'),
		out_completeness = os.path.join(evaluation_output_folder, 'merqury.completeness.stats')
	params:
		merqury_prefix = os.path.join(evaluation_output_folder, "merqury"),
		asm_a = os.path.join(config['results'], prefix, "contigging/hifiasm", prefix + '.asm.a_ctg.fasta') if (contig_assembler == "hifiasm" and not phasing_mode) else "" 
	conda:
		os.path.join(config['snakemake_dir_path'], 'envs/evaluation.yml')
	threads:
		resource['merqury']['threads']
	resources:
		mem_mb=resource['merqury']['mem_mb'],
		time=resource['merqury']['time']
	log:
		os.path.join(config['results'], "logs", prefix, 'merqury.log')
	shell:
		"(merqury.sh {input.meryl_db} {input.asm_p} {params.asm_a} {params.merqury_prefix}) 2> {log};"

rule busco:
	""" Execute BUSCO following Hifiasm to assess the completeness of the assembled genome based on conserved genes."""
		
	input:
		asm_p = contig_output_primary(contig_assembler, phasing_mode)
	output:
		outdir = directory(os.path.join(evaluation_output_folder, 'busco')),
		busco_out = os.path.join(evaluation_output_folder, 'busco', 'short_summary.specific.' + config['busco_lineage'] + '_odb10.busco_out.txt')
	params:
		out_path = evaluation_output_folder,
		out_name = 'busco',
		lineage = config['busco_lineage']
	conda:
		os.path.join(config['snakemake_dir_path'], 'envs/evaluation.yml')
	threads:
		resource['busco']['threads']
	resources:
		mem_mb=resource['busco']['mem_mb'],
		time=resource['busco']['time']
	log:
		os.path.join(config['results'], "logs", prefix, 'busco.log')
	shell:
		"(busco -f -m genome -i {input.asm_p} -o {params.out_name} --out_path {params.out_path} -l {params.lineage} -c {threads}) 2> {log};"


rule gfastats:
	""" run gfastats on assembly"""

	input:
		asm_p = contig_output_primary(contig_assembler, phasing_mode)
	output:
		out_gfastats = os.path.join(evaluation_output_folder, "stats", 'gfastats.txt')
	conda:
		os.path.join(config['snakemake_dir_path'], 'envs/evaluation.yml')
	threads:
		resource['gfastats']['threads']
	resources:
		mem_mb=resource['gfastats']['mem_mb'],
		time=resource['gfastats']['time']
	log:
		os.path.join(config['results'], "logs", prefix, 'gfastats.log')
	shell:
		"(gfastats --nstar-report {input.asm_p} > {output}) 2> {log};"


rule merqury_purging:
	""" Execute Merqury following Hifiasm to evaluate the quality and accuracy of the assembled genome."""
	input:
		asm_p = os.path.join(config['results'], prefix, "contigging/hifiasm_hic_l{l}", prefix + '.asm.hic.p_ctg.fa'),
		meryl_db = os.path.join(config['results'] , prefix, "genome_profiling", (prefix +".meryl"))
	output:
		out_qv = os.path.join(config['results'], prefix, "assembly_evaluation/hifiasm_hic_l{l}/merqury.qv"),
		out_completeness = os.path.join(config['results'], prefix, "assembly_evaluation/hifiasm_hic_l{l}/merqury.completeness.stats")
	params:
		merqury_prefix = os.path.join(config['results'], prefix, "assembly_evaluation/hifiasm_hic_l{l}/merqury"),
		asm_a = os.path.join(config['results'], prefix, "contigging/hifiasm", prefix + '.asm.a_ctg.fasta') if (contig_assembler == "hifiasm" and not phasing_mode) else "" 
	conda:
		os.path.join(config['snakemake_dir_path'], 'envs/evaluation.yml')
	threads:
		resource['merqury']['threads']
	resources:
		mem_mb=resource['merqury']['mem_mb'],
		time=resource['merqury']['time']
	log:
		os.path.join(config['results'], "logs", prefix, 'merqury_hifiasm_l{l}.log')
	shell:
		"(merqury.sh {input.meryl_db} {input.asm_p} {params.asm_a} {params.merqury_prefix}) 2> {log};"


rule busco_purging:
	""" Execute BUSCO following Hifiasm to assess the completeness of the assembled genome based on conserved genes."""
		
	input:
		asm_p = os.path.join(config['results'], prefix, "contigging/hifiasm_hic_l{l}", prefix + '.asm.hic.p_ctg.fa'),
	output:
		outdir = directory(os.path.join(config['results'], prefix, "assembly_evaluation/hifiasm_hic_l{l}/busco"),),
		busco_out = os.path.join(config['results'], prefix, "assembly_evaluation/hifiasm_hic_l{l}/busco", 'short_summary.specific.' + config['busco_lineage'] + '_odb10.busco_out.txt')
	params:
		out_path = os.path.join(config['results'], prefix, "assembly_evaluation/hifiasm_hic_l{l}"),
		out_name = 'busco',
		lineage = config['busco_lineage']
	conda:
		os.path.join(config['snakemake_dir_path'], 'envs/evaluation.yml')
	threads:
		resource['busco']['threads']
	resources:
		mem_mb=resource['busco']['mem_mb'],
		time=resource['busco']['time']
	log:
		os.path.join(config['results'], "logs", prefix, 'busco_purging_l{l}.log')
	shell:
		"(busco -f -m genome -i {input.asm_p} -o {params.out_name} --out_path {params.out_path} -l {params.lineage} -c {threads}) 2> {log};"


rule gfastats_purging:
	""" run gfastats on assembly"""

	input:
		asm_p = os.path.join(config['results'], prefix, "contigging/hifiasm_hic_l{l}", prefix + '.asm.hic.p_ctg.fa'),
	output:
		out_gfastats = os.path.join(config['results'], prefix, "assembly_evaluation/hifiasm_hic_l{l}", "stats", 'gfastats.txt')
	conda:
		os.path.join(config['snakemake_dir_path'], 'envs/evaluation.yml')
	threads:
		resource['gfastats']['threads']
	resources:
		mem_mb=resource['gfastats']['mem_mb'],
		time=resource['gfastats']['time']
	log:
		os.path.join(config['results'], "logs", prefix, 'gfastats_purging_l{l}.log')
	shell:
		"(gfastats --nstar-report {input.asm_p} > {output}) 2> {log};"

