configfile: 'seekdeep_illumina_general.yaml'

def get_qluster_fastqs(wildcards):
	output_files=[]
	samples=[line.strip() for line in open(config['output_folder']+'/non-empty_extractions.txt')]
	for sample in samples:
		sample_prefix=sample.split('MID')[1].split('_')[0]
		output_files.append(config['output_folder']+f'/analysis/{sample_prefix}_extraction/{sample}_qlusterOut/output.fastq.gz')
	return output_files

import os
amplicon_folder=config['output_folder']+'/analysis/popClustering'
amplicons=os.listdir(amplicon_folder)
amplicons.remove('locationByIndex')

rule all:
	input:
#		all_sample_commands=expand(config['output_folder']+'/qluster_shell_commands/{sample}_qluster_command.sh', sample=[line.strip() for line in open(config['output_folder']+'/non-empty_extractions.txt')]),
#		qluster_output=get_qluster_fastqs
#		all_process_cluster_commands=expand(config['output_folder']+'/process_cluster_shell_commands/{amplicon}_process_cluster_command.sh', amplicon=amplicons)
		analysis_folder=expand(config['output_folder']+'/analysis/popClustering/{amplicon}/analysis/selectedClustersInfo.tab.txt.gz', amplicon=amplicons)


rule prep_qluster:
	input:
		runnable_samples=config['output_folder']+'/non-empty_extractions.txt',
		qluster_commands=config['output_folder']+'/analysis/qlusterCmds.txt'
	params:
		output_folder=config['output_folder']+'/qluster_shell_commands'
	output:
		all_sample_commands=expand(config['output_folder']+'/qluster_shell_commands/{sample}_qluster_command.sh', sample=[line.strip() for line in open(config['output_folder']+'/non-empty_extractions.txt')])
	script:
		'scripts/prep_qluster.py'

rule run_qluster:
	input:
		sample_command=config['output_folder']+'/qluster_shell_commands/{sample}_qluster_command.sh',
		data_folder=config['primer_plus_fastq_binding'],
		sif_file=config['sif_file_location'],
		genome_root_folder=config['genome_binding'],
		actual_qluster_file=config['output_folder']+'/qluster_shell_commands/{sample}_qluster_command.sh'
	params:
		output_dir=config['output_folder'],
		softlink_fastq_binding=config['softlink_fastq_binding'],
		singularity_qluster_file='/seekdeep_output/qluster_shell_commands/{sample}_qluster_command.sh'
	output:
		qluster_output=config['output_folder']+'/analysis/{sample_prefix}_extraction/{sample}_qlusterOut/output.fastq.gz'
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
		{input.sif_file} bash {params.singularity_qluster_file}
		'''

rule prep_process_cluster:
	input:
		process_cluster_commands=config['output_folder']+'/analysis/processClusterCmds.txt',
		qluster_done=get_qluster_fastqs
	params:
		output_folder=config['output_folder']+'/process_cluster_shell_commands'
	output:
		all_process_cluster_commands=expand(config['output_folder']+'/process_cluster_shell_commands/{amplicon}_process_cluster_command.sh', amplicon=amplicons)
	script:
		'scripts/prep_process_cluster.py'

rule run_process_cluster:
	input:
		data_folder=config['primer_plus_fastq_binding'],
		sif_file=config['sif_file_location'],
		genome_root_folder=config['genome_binding'],
		actual_process_cluster_command=config['output_folder']+'/process_cluster_shell_commands/{amplicon}_process_cluster_command.sh'
	params:
		output_dir=config['output_folder'],
		softlink_fastq_binding=config['softlink_fastq_binding'],
		singularity_process_cluster_command='/seekdeep_output/process_cluster_shell_commands/{amplicon}_process_cluster_command.sh'
	output:
		analysis_folder=config['output_folder']+'/analysis/popClustering/{amplicon}/analysis/selectedClustersInfo.tab.txt.gz'
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
		{input.sif_file} bash {params.singularity_process_cluster_command}
		'''

#gen_config_commands=[line.strip() for line in open(config['output_folder']+'/analysis/genConfigCmds.txt')]
#rule gen_config:
#	input:
#		data_folder=config['primer_plus_fastq_binding'],
#		sif_file=config['sif_file_location'],
#		setup_done=config['output_folder']+'/finished_setup.txt',
#		genome_root_folder=config['genome_binding'],
#		process_cluster_files=expand(config['output_folder']+'/process_cluster_jobs/{number}_process_cluster_done.txt', number=list(range(len(process_cluster_commands))))
#	params:
#		output_dir=config['output_folder'],
#		softlink_fastq_binding=config['softlink_fastq_binding'],
#		command=lambda wildcards: gen_config_commands[int(wildcards.number)]
#	output:
#		gen_config_done=config['output_folder']+'/gen_config_jobs/{number}_gen_config_done.txt',
#	threads: config['cpus_to_use']
#	resources:
#		time_min=config['max_run_time_min'],
#		mem_mb=config['max_memory_mb'],
#		nodes=config['cpus_to_use']
#	shell:
#		'''
#		singularity exec -B {input.data_folder}:/input_data \
#		-B {params.output_dir}:/seekdeep_output \
#		-B {input.genome_root_folder}:/genome_info \
#		{params.softlink_fastq_binding} \
#		-H {params.output_dir}/analysis/:/home/analysis \
#		{input.sif_file} {params.command}
#		touch {output.gen_config_done}
#		'''

#rule all_done:
#	input:
#		gen_config_files=expand(config['output_folder']+'/gen_config_jobs/{number}_gen_config_done.txt', number=list(range(len(gen_config_commands))))
#	output:
#		all_finished=config['output_folder']+'/all_finished.txt'
#	shell:
#		'''
#		touch {output.all_finished}
#		'''
