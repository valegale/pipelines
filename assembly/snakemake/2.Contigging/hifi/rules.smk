rule hifiasm:
	"""Execute Hifiasm to generate primary and alternative assemblies."""

	input:
		hifi_fastq = expand("{sample}", sample=samples['file_path'])
	output:
		asm_p_fa = os.path.join(config['results'], '2.Contigging/hifiasm/purgeL' + config['purgel'] + '/hifiasmL' + config['purgel'] + '.asm.p_ctg.gfa'),
		asm_a_fa = os.path.join(config['results'], '2.Contigging/hifiasm/purgeL' + config['purgel'] + '/hifiasmL' + config['purgel'] + '.asm.a_ctg.gfa')
	log:
		os.path.join(config['snakemake_dir_path'] + 'logs/hifiasmL' + config['purgel'] + '.log')
	conda:
		os.path.join(config['snakemake_dir_path'], 'envs/hifiasm.yaml')
	threads:
		16
	resources:
		mem_mb = 70000
	params:
		runtime = '40:00:00',
		prefix = os.path.join(config['results'], '2.Contigging/hifiasm/purgeL' + config['purgel'] + '/hifiasmL' + config['purgel'] + '.asm'),
		purgel=config['purgel']
	shell:
		"""
		(hifiasm -o {params.prefix}  --primary -t {threads} -l {params.purgel} {input.hifi_fastq}) 2> {log}
		"""


rule hicanu:
	"""Execute HiCanu to generate primary and alternative assemblies."""
	input:
		hifi_fastq = expand("{sample}", sample=samples['file_path'])
	output:
		outdir = directory(os.path.join(config['results'], '2.Contigging/hicanu')),
		outasm = os.path.join(config['results'], '2.Contigging/hicanu/asm.contigs.fasta')
	log:
		os.path.join(config['snakemake_dir_path'], 'logs/hicanu.log')
	conda:
		os.path.join(config['snakemake_dir_path'], 'envs/canu.yaml')
	threads:
		16
	resources:
		mem_gb = 70
	params:
		runtime = '40:00:00',
		genomeSize = config['genome_size_estimate']
	shell:
		"""
		(canu -p asm -d {output.outdir} genomeSize={params.genomeSize} -pacbio-hifi {input.hifi_fastq} useGrid=false Threads {threads} Memory {resources.mem_gb}) 2> {log}
		"""

rule flye:
	"""Execute Flye to generate primary and alternative assemblies."""

	input:
		hifi_fastq = expand("{sample}", sample=samples['file_path'])
	output:
		outdir = directory(os.path.join(config['results'], '2.Contigging/flye')),
		asm_flye = os.path.join(config['results'], '2.Contigging/flye/assembly.fasta')
	log:
		os.path.join(config['snakemake_dir_path'], 'logs/2_Build_asm/flye/flye.log')
	threads:
		16
	resources:
		mem_mb = 70000
	conda:
		os.path.join(config['snakemake_dir_path'], "envs/flye.yaml")
	params:
		genomeSize = config['genome_size_estimate']
	shell:
		"""
		(flye --threads {threads} --out-dir {output.outdir} --genome-size {params.genomeSize} --pacbio-hifi {input.hifi_fastq}) 2> {log}
		"""

rule gfa_to_fasta:
	"""Convert the GFA output from Hifiasm into a FASTA format."""
	input:
		asm_p_gfa = os.path.join(config['results'], '2.Contigging/hifiasm/purgeL' + config['purgel'] + '/hifiasmL' + config['purgel'] + '.asm.p_ctg.gfa'),
		asm_a_gfa = os.path.join(config['results'], '2.Contigging/hifiasm/purgeL' + config['purgel'] + '/hifiasmL' + config['purgel'] + '.asm.a_ctg.gfa'),
	output:
		asm_p_fa = os.path.join(config['results'], '2.Contigging/hifiasm/purgeL'  + config['purgel'] + '/hifiasmL' + config['purgel'] + '.asm.p_ctg.fasta'),
		asm_a_fa = os.path.join(config['results'], '2.Contigging/hifiasm/purgeL'  + config['purgel'] + '/hifiasmL' + config['purgel'] + '.asm.a_ctg.fasta')
	log:
		os.path.join(config['snakemake_dir_path'], 'logs/gfa_to_fa_L',config['purgel'],'.log')
	threads:
		16
	resources:
		mem_mb = 1000
	params:
		ploidy = config['ploidy'],
		runtime = '03:00:00'
	run:
		shell("""awk '/^S/{{print ">"$2;print $3}}' {input.asm_p_gfa} > {output.asm_p_fa} 2> {log};
			awk '/^S/{{print ">"$2;print $3}}' {input.asm_a_gfa} > {output.asm_a_fa} 2>> {log}
			""")
