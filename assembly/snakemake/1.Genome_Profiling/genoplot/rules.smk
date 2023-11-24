rule count_kmers:
	input:
		reads=lambda wildcards: samples[samples["sample_nr"] == wildcards.sample_nr]["file_path"]
	output:
		directory(os.path.join(config['results'],"1.Genome_Profiling/{sample_nr}.meryl"))
	conda:
		os.path.join(config['snakemake_dir_path'], "envs/genome_profiling.yaml")
	threads:
		3
	shell:
		"meryl count k={k} memory=10G threads={threads} {input.reads} output {output}"

rule merge_kmers:
	input:
		expand(os.path.join(config['results'],"1.Genome_Profiling/{sample_nr}.meryl"), sample_nr=samples['sample_nr']),  
	output:
		directory(output_meryl)
	conda:
		os.path.join(config['snakemake_dir_path'], "envs/genome_profiling.yaml")
	threads:
		6
	shell:
		"meryl union-sum memory=10G threads={threads} output {output} {input}"

rule create_hist:
	input:
		output_meryl
	output:
		os.path.join(config['results'],"1.Genome_Profiling", (config["prefix"]+".hist"))
	conda:
		os.path.join(config['snakemake_dir_path'], "envs/genome_profiling.yaml")
	threads:
		1
	shell:
		"meryl histogram {input} > {output}"

rule run_genomescope:
	input:
		os.path.join(config['results'],"1.Genome_Profiling", (config["prefix"]+".hist"))
	output:
		directory(os.path.join(config['results'],"1.Genome_Profiling", (config["prefix"]+"_genomescope")))
	conda:
		os.path.join(config['snakemake_dir_path'], "envs/genome_profiling.yaml")
	threads:
		1
	shell:
		"genomescope2 -i {input} -p {ploidy} -k {k} -o {output}"

rule get_genomescopeStats:
	input:
		os.path.join(config['results'],"1.Genome_Profiling", (config["prefix"]+"_genomescope"))
	output:
		"Estimated_genome_size"
	threads:
		1
	shell: # change here
		"bash 1.Genome_Profiling/genoplot/get_stats.sh {input}"

rule prepare_kmersSmudge:
	input:
		output_meryl
	output:
		os.path.join(config['results'],"1.Genome_Profiling", (str(prefix)+"_L"+str(L)+"_U"+str(U)+".kmers"))
	conda:
		os.path.join(config['snakemake_dir_path'], "envs/genome_profiling.yaml")
	threads:
		6
	shell:
		"meryl print less-than 1500 greater-than 10 threads={threads} memory=10G {input} | sort > {output}"

rule run_smudgeplot:
	input:
		os.path.join(config['results'],"1.Genome_Profiling", (str(prefix)+"_L"+str(L)+"_U"+str(U)+".kmers"))
	output:
		os.path.join(config['results'],"1.Genome_Profiling", (str(prefix)+"_L"+str(L)+"_U"+str(U)+"_coverages.tsv"))
	conda:
		os.path.join(config['snakemake_dir_path'], "envs/genome_profiling.yaml")
	threads:
		1
	shell:
		"smudgeplot.py hetkmers -o {prefix}_L{L}_U{U} < {input}"

rule create_smudgeplot:
	input:
		os.path.join(config['results'],"1.Genome_Profiling", (str(prefix)+"_L"+str(L)+"_U"+str(U)+"_coverages.tsv"))
	output:
		output_smudgeplot
	conda:
		os.path.join(config['snakemake_dir_path'], "envs/genome_profiling.yaml")
	threads:
		1
	shell:
		"smudgeplot_plot.R -L {L} -i {input} -o {prefix}_L{L}_U{U}"