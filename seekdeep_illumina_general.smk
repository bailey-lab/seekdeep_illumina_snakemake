configfile: 'seekdeep_illumina_general.yaml'
rule all:
	input:
		analysis_done=config['output_folder']+'/finished_analysis.txt',
		out_snakefile=config['output_folder']+'/seekdeep_illumina_general.smk',
		out_config_file=config['output_folder']+'/seekdeep_illumina_general.yaml'

rule copy_files:
	'''
	copies snakemake script and config files to output folder for reproducibility.
	'''
	input:
		in_snakefile='seekdeep_illumina_general.smk',
		in_config_file='seekdeep_illumina_general.yaml'
	output:
		out_snakefile=config['output_folder']+'/seekdeep_illumina_general.smk',
		out_config_file=config['output_folder']+'/seekdeep_illumina_general.yaml'
	shell:
		'''
		cp {input.in_snakefile} {output.out_snakefile}
		cp {input.in_config_file} {output.out_config_file}
		'''

rule genTargetInfoFromGenomes:
	input:
		data_folder=config['primer_plus_fastq_binding'],
		genome_root_folder=config['genome_binding'],
		sif_file=config['sif_file_location']
	params:
		output_dir=config['output_folder'],
		primer_file=config['primer_file'],
		gff_subfolder=config['gff_subfolder'],
		genome_subfolder=config['genome_subfolder'],
		read_length=config['read_length'],
		extra_args=config['extra_gen_target_info_cmds']
	output:
		primer_info=config['output_folder']+'/extractedRefSeqs/locationsByGenome/Pf3D7_infos.tab.txt'
	threads: config['cpus_to_use']
	resources:
		time_min=config['max_run_time_min'],
		mem_mb=config['max_memory_mb'],
		nodes=config['cpus_to_use']
	shell:
		'''
		singularity exec -B {input.genome_root_folder}:/genome_info \
		-B {input.data_folder}:/input_data \
		-B {params.output_dir}:/seekdeep_output {input.sif_file} \
		SeekDeep genTargetInfoFromGenomes --primers /input_data/{params.primer_file} \
		--pairedEndLength {params.read_length} --genomeDir /genome_info/{params.genome_subfolder} \
		--gffDir /genome_info/{params.gff_subfolder} {params.extra_args} \
		--dout /seekdeep_output/extractedRefSeqs --overWriteDir --numThreads \
		{threads}
		'''

rule setupTarAmpAnalysis:
	input:
		data_folder=config['primer_plus_fastq_binding'],
		genome_root_folder=config['genome_binding'],
		sif_file=config['sif_file_location'],
		primer_info=config['output_folder']+'/extractedRefSeqs/locationsByGenome/Pf3D7_infos.tab.txt'
	params:
		output_dir=config['output_folder'],
		primer_file=config['primer_file'],
		fastq_folder=config['fastq_subfolder'],
		genome_subfolder=config['genome_subfolder'],
		sample_names=config['sample_names'],
		for_seekdeep='/seekdeep_output/extractedRefSeqs/forSeekDeep',
		softlink_fastq_binding=config['softlink_fastq_binding'],
		extra_args=config['extra_setup_tar_amp_cmds'],
		extra_extractor_cmds=config['extra_extractor_cmds'],
		extra_qluster_cmds=config['extra_qluster_cmds'],
		extra_process_cluster_cmds=config['extra_process_cluster_cmds']
	threads: config['cpus_to_use']
	output:
		setup_done=config['output_folder']+'/finished_setup.txt'
	shell:
		'''
		singularity exec -B {input.data_folder}:/input_data \
		-B {input.genome_root_folder}:/genome_info \
		-B {params.output_dir}:/seekdeep_output \
		{params.softlink_fastq_binding} {input.sif_file} \
		SeekDeep setupTarAmpAnalysis --samples /input_data/{params.sample_names} \
		--outDir /seekdeep_output/analysis \
		--inputDir /input_data/{params.fastq_folder} \
		--idFile /input_data/{params.primer_file} --lenCutOffs \
		{params.for_seekdeep}/lenCutOffs.txt \
		--overlapStatusFnp {params.for_seekdeep}/overlapStatuses.txt \
		--refSeqsDir {params.for_seekdeep}/refSeqs/ {params.extra_args} \
		{params.extra_extractor_cmds} {params.extra_qluster_cmds} \
		{params.extra_process_cluster_cmds} --numThreads {threads}
		touch {output.setup_done}
		'''

extractor_commands=[line.strip() for line in open(config['output_folder']+'/analysis/extractorCmds.txt')]
rule run_extractor:
	input:
		data_folder=config['primer_plus_fastq_binding'],
		sif_file=config['sif_file_location'],
		setup_done=config['output_folder']+'/finished_setup.txt',
		genome_root_folder=config['genome_binding']
	params:
		output_dir=config['output_folder'],
		softlink_fastq_binding=config['softlink_fastq_binding'],
		command=lambda wildcards: extractor_commands[int(wildcards.number)]
	output:
		extractor_done=config['output_folder']+'/extractor_jobs/{number}_extractor_done.txt'
	threads: config['cpus_to_use']
	resources:
		time_min=config['max_run_time_min'],
		mem_mb=config['max_memory_mb'],
		nodes=config['cpus_to_use']
	shell:
		'''
		singularity exec -B {input.data_folder}:/input_data \
		-B {params.output_dir}:/seekdeep_output \
		-B {input.genome_root_folder}:/genome_info \
		{params.softlink_fastq_binding} \
		-H {params.output_dir}/analysis/:/home/analysis \
		{input.sif_file} {params.command}
		touch {output.extractor_done}
		'''

qluster_commands=[line.strip() for line in open(config['output_folder']+'/analysis/qlusterCmds.txt')]
rule run_qluster:
	input:
		data_folder=config['primer_plus_fastq_binding'],
		sif_file=config['sif_file_location'],
		genome_root_folder=config['genome_binding']
		extractor_files=expand(config['output_folder']+'/extractor_jobs/{number}_extractor_done.txt', number=list(range(extractor_commands)))
	params:
		output_dir=config['output_folder'],
		softlink_fastq_binding=config['softlink_fastq_binding'],
		command=lambda wildcards: qluster_commands[int(wildcards.number)]
	output:
		qluster_done=config['output_folder']+'/qluster_jobs/{number}_qluster_done.txt')
	threads: config['cpus_to_use']
	resources:
		time_min=config['max_run_time_min'],
		mem_mb=config['max_memory_mb'],
		nodes=config['cpus_to_use']
	shell:
		'''
		singularity exec -B {input.data_folder}:/input_data \
		-B {params.output_dir}:/seekdeep_output \
		-B {input.genome_root_folder}:/genome_info \
		{params.softlink_fastq_binding} \
		-H {params.output_dir}/analysis/:/home/analysis \
		{input.sif_file} {params.command}
		touch {output.qluster_done}
		'''

process_cluster_commands=[line.strip() for line in open(config['output_folder']+'/analysis/processClusterCmds.txt')]
rule run_process_cluster:
	input:
		data_folder=config['primer_plus_fastq_binding'],
		sif_file=config['sif_file_location'],
		setup_done=config['output_folder']+'/finished_setup.txt',
		genome_root_folder=config['genome_binding']
		qluster_files=expand(config['output_folder']+'/qluster_jobs/{number}_qluster_done.txt', number=list(range(qluster_commands)))
	params:
		output_dir=config['output_folder'],
		softlink_fastq_binding=config['softlink_fastq_binding'],
		command=lambda wildcards: process_cluster_commands[int(wildcards.number)]
	output:
		process_cluster_done=config['output_folder']+'/process_cluster_jobs/{number}_process_cluster_done.txt'
	threads: config['cpus_to_use']
	resources:
		time_min=config['max_run_time_min'],
		mem_mb=config['max_memory_mb'],
		nodes=config['cpus_to_use']
	shell:
		'''
		singularity exec -B {input.data_folder}:/input_data \
		-B {params.output_dir}:/seekdeep_output \
		-B {input.genome_root_folder}:/genome_info \
		{params.softlink_fastq_binding} \
		-H {params.output_dir}/analysis/:/home/analysis \
		{input.sif_file} {params.command}
		touch {output.process_cluster_done}
		'''

gen_config_commands=[line.strip() for line in open(config['output_folder']+'/analysis/genConfigCmds.txt')]
rule gen_config:
	input:
		data_folder=config['primer_plus_fastq_binding'],
		sif_file=config['sif_file_location'],
		setup_done=config['output_folder']+'/finished_setup.txt',
		genome_root_folder=config['genome_binding']
		process_cluster_files=expand(config['output_folder']+'/process_cluster_jobs/{number}_process_cluster_done.txt', number=list(range(process_cluster_commands)))
	params:
		output_dir=config['output_folder'],
		softlink_fastq_binding=config['softlink_fastq_binding'],
		command=lambda wildcards: gen_config_commands[int(wildcards.number)]
	output:
		gen_config_done=config['output_folder']+'/gen_config_jobs/{number}_gen_config_done.txt'
	threads: config['cpus_to_use']
	resources:
		time_min=config['max_run_time_min'],
		mem_mb=config['max_memory_mb'],
		nodes=config['cpus_to_use']
	shell:
		'''
		singularity exec -B {input.data_folder}:/input_data \
		-B {params.output_dir}:/seekdeep_output \
		-B {input.genome_root_folder}:/genome_info \
		{params.softlink_fastq_binding} \
		-H {params.output_dir}/analysis/:/home/analysis \
		{input.sif_file} {params.command}
		touch {output.gen_config_done}
		'''

