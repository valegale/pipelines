rule count_kmers:
	input:
		reads=lambda wildcards: samples[samples["sample_nr"] == wildcards.sample_nr]["file_path"]
	output:
		directory("{sample_nr}.meryl")
	conda:
		"genoplots"
	threads:
		3
	shell:
		"meryl count k={k} memory=10G threads={threads} {input.reads} output {output}"

rule merge_kmers:
	input:
		expand("{sample_nr}.meryl", sample_nr=samples['sample_nr']),  
	output:
		directory(output_meryl)
	conda:
		"genoplots"
	threads:
		6
	shell:
		"meryl union-sum memory=10G threads={threads} output {output} {input}"

rule create_hist:
	input:
		output_meryl
	output:
		config["prefix"]+".hist"
	conda:
		"genoplots"
	threads:
		1
	shell:
		"meryl histogram {input} > {output}"

rule run_genomescope:
	input:
		config["prefix"]+".hist"
	output:
		directory(config["prefix"]+"_genomescope")
	conda:
		"genoplots"
	threads:
		1
	shell:
		"genomescope2 -i {input} -p {ploidy} -k {k} -o {output}"

rule get_genomescopeStats:
	input:
		config["prefix"]+"_genomescope"
	output:
		"Estimated_genome_size"
	threads:
		1
	shell:
		"bash get_stats.sh {input}"

rule prepare_kmersSmudge:
	input:
		output_meryl
	output:
		str(prefix)+"_L"+str(L)+"_U"+str(U)+".kmers"
	conda:
		"genoplots"
	threads:
		6
	shell:
		"meryl print less-than 1500 greater-than 10 threads={threads} memory=10G {input} | sort > {output}"

rule run_smudgeplot:
	input:
		str(prefix)+"_L"+str(L)+"_U"+str(U)+".kmers"
	output:
		str(prefix)+"_L"+str(L)+"_U"+str(U)+"_coverages.tsv"
	conda:
		"genoplots"
	threads:
		1
	shell:
		"smudgeplot.py hetkmers -o {prefix}_L{L}_U{U} < {input}"

rule create_smudgeplot:
	input:
		str(prefix)+"_L"+str(L)+"_U"+str(U)+"_coverages.tsv"
	output:
		output_smudgeplot
	conda:
		"genoplots"
	threads:
		1
	shell:
		"smudgeplot_plot.R -L {L} -i {input} -o {prefix}_L{L}_U{U}"