configfile: 'seekdeep_illumina_general.yaml'

rule all:
	input:
#		profile=expand(config['output_folder']+'/analysis/{sample}_extraction/extractionProfile.tab.txt', sample=config['samples'])
		runnable_samples=config['output_folder']+'/non-empty_extractions.txt',

rule prep_extractor:
	input:
		extractor_commands=config['output_folder']+'/analysis/extractorCmds.txt'
	params:
		output_folder=config['output_folder']+'/extractor_shell_commands'
	output:
		all_sample_commands=expand(config['output_folder']+'/extractor_shell_commands/{sample}_extraction_command.sh', sample=config['samples'])
	script:
		'scripts/prep_extractor.py'

rule run_extractor:
	input:
		data_folder=config['primer_plus_fastq_binding'],
		sif_file=config['sif_file_location'],
#		setup_done=config['output_folder']+'/finished_setup.txt',
		genome_root_folder=config['genome_binding'],
		actual_shell_script=config['output_folder']+'/extractor_shell_commands/{sample}_extraction_command.sh'	
	params:
		output_dir=config['output_folder'],
		softlink_fastq_binding=config['softlink_fastq_binding'],
		singularity_shell_script='/seekdeep_output/extractor_shell_commands/{sample}_extraction_command.sh'	
	output:
		profile=config['output_folder']+'/analysis/{sample}_extraction/extractionProfile.tab.txt',
		folder=directory(config['output_folder']+'/analysis/{sample}_extraction')
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
		{input.sif_file} bash {params.singularity_shell_script}
		'''

rule analyze_extractor:
	input:
		profiles=expand(config['output_folder']+'/analysis/{sample}_extraction/extractionProfile.tab.txt', sample=config['samples']),
		folders=expand(config['output_folder']+'/analysis/{sample}_extraction', sample=config['samples'])
	output:
		runnable_samples=config['output_folder']+'/non-empty_extractions.txt'
	script:
		'scripts/analyze_extractor.py'
