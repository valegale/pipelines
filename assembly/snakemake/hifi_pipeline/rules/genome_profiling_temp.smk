### meryl

rule count_kmers:
	input:
		reads=lambda wildcards: samples[samples["sample_nr"] == wildcards.sample_nr]["hifi_data"]
	output:
		temp(directory(os.path.join(config['results'],"genome_profiling/{sample_nr}.meryl")))
	conda:
		os.path.join(config['snakemake_dir_path'], "envs/genome_profiling.yml")
	threads:
		resource['meryl_count']['threads']
	resources:
		mem_mb=resource['meryl_count']['mem_mb'],
		time=resource['meryl_count']['time']
	shell:
		"meryl count k={kmer} memory={resources.mem_mb} threads={threads} {input.reads} output {output}"

rule merge_kmers:
	input:
		expand(os.path.join(config['results'],"genome_profiling/{sample_nr}.meryl"), sample_nr=samples['sample_nr']),  
	output:
		directory(os.path.join(config['results'],"genome_profiling", (config["prefix"] +".meryl")))
	conda:
		os.path.join(config['snakemake_dir_path'], "envs/genome_profiling.yml")
	threads:
		resource['meryl_merge']['threads']
	resources:
		mem_mb=resource['meryl_merge']['mem_mb'],
		time=resource['meryl_merge']['time']
	shell:
		"meryl union-sum memory={resources.mem_mb} threads={threads} output {output} {input} "

rule create_hist:
	input:
		os.path.join(config['results'],"genome_profiling", (config["prefix"] +".meryl"))
	output:
		os.path.join(config['results'],"genome_profiling", (prefix + ".hist"))
	conda:
		os.path.join(config['snakemake_dir_path'], "envs/genome_profiling.yml")
	threads:
		resource['meryl_histogram']['threads']
	resources:
		mem_mb=resource['meryl_histogram']['mem_mb'],
		time=resource['meryl_histogram']['time']
	shell:
		"meryl histogram {input} > {output}"

#### smudgeplot

rule run_smudgeplot:
	input:
		hist=os.path.join(config['results'],"genome_profiling", (prefix + ".hist"))
	output:
		plot1=os.path.join(config['results'],"genome_profiling/{prefix}/smudgeplot/results_smudgeplot.png"),
		plot2=os.path.join(config['results'],"genome_profiling/{prefix}/smudgeplot/results_smudgeplot_log10.png"),
		summary_table=os.path.join(config['results'],"genome_profiling/{prefix}/smudgeplot/results_summary_table.tsv")
	params:
		output_path=os.path.join(config['results'],"genome_profiling/{prefix}/smudgeplot"),
		meryl_path=os.path.join(config['results'],"genome_profiling", (config["prefix"] +".meryl")),
		kmer=kmer
	conda:
		os.path.join(config['snakemake_dir_path'], "envs/genome_profiling.yml")
	threads:
		resource['run_smudgeplot']['threads']
	resources:
		mem_mb=resource['run_smudgeplot']['mem_mb'],
		time=resource['run_smudgeplot']['time']
	log:
		os.path.join(config['results'],"logs/{prefix}_smudgeplot.log")
	shell:
		"""
		L=$(smudgeplot.py cutoff {input.hist} L)
		U=$(smudgeplot.py cutoff {input.hist} U)

		(meryl print less-than ${{U}} greater-than ${{L}} threads={threads} memory={resources.mem_mb} {params.meryl_path} | sort > {params.output_path}/meryl_L${{L}}_U${{U}}.dump) &> {log}

		(smudgeplot.py hetkmers -o {params.output_path}/meryl_L${{L}}_U${{U}} --middle {params.output_path}/meryl_L${{L}}_U${{U}}.dump) &> {log}

		(smudgeplot.py plot -o {params.output_path}/results {params.output_path}/meryl_L${{L}}_U${{U}}_coverages.tsv -k {params.kmer} ) &> {log}
		"""

#### genomescope, still testing!

# check if prefix and ploidy are taken automatically
rule run_genomescope:
	input:
		os.path.join(config['results'],"genome_profiling", (prefix + ".hist"))
	output:
    	directory(os.path.join(config['results'],"genome_profiling/{prefix}/genomescope"))
	conda:
		os.path.join(config['snakemake_dir_path'], "envs/genome_profiling.yml")
	threads:
		resource['run_genomescope']['threads']
	resources:
		mem_mb=resource['run_genomescope']['mem_mb'],
		time=resource['run_genomescope']['time']
	shell:
		"genomescope2 -i {input} -p {ploidy} -k {k} -o {output}"

# change the script, creating two files or just one with statistics? 
rule get_genomescope_stats:
	input:
		os.path.join(config['results'],"genome_profiling/{prefix}/genomescope")
	output:
		"Estimated_genome_size"
	threads:
		1
	shell: # change here
		"bash 1.Genome_Profiling/genoplot/get_stats.sh {input}"
