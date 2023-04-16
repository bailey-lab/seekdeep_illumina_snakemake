configfile: 'seekdeep_nanopore_general.yaml'

rule all:
	input:
#		primer_info=config['output_folder']+'/extractedRefSeqs/locationsByGenome/Pf3D7_infos.tab.txt'
#		done=config['output_folder']+'/finished_setup.txt'
		analysis_done=config['output_folder']+'/finished_analysis.txt'	
#		analysis_done=config['output_folder']+'/analysis/popClustering'

rule copy_files:
	'''
	copies snakemake script and config files to output folder for reproducibility.
	'''
	input:
		in_snakefile='seekdeep_illumina_general.smk'
		in_config_file='seekdeep_illumina_general.yaml'
	output:
		out_snakefile=config['output_folder']+'seekdeep_nanopore_general.smk',
		out_config_file=config['output_folder']+'seekdeep_nanopore_general.yaml'
	shell:
		'''
		cp {input.in_snakefile} {output.out_snakefile}
		cp {input.in_config_file} {output.out_config_file}
		'''

rule get_primer_info:
	input:
		data_folder=config['primer_plus_fastq_binding'],
		genome_root_folder=config['genome_binding'],
		sif_file=config['sif_file_location']
	params:
		output_dir=config['output_folder'],
		primer_file=config['primer_file'],
		gff_subfolder=config['gff_subfolder'],
		genome_subfolder=config['genome_subfolder'],
		insert_size=config['insert_size']
	output:
		primer_info=config['output_folder']+'/extractedRefSeqs/locationsByGenome/Pf3D7_infos.tab.txt'
	threads: config['cpus_to_use']
	shell:
		'''
		singularity exec -B {input.genome_root_folder}:/genome_info \
		-B {input.data_folder}:/input_data \
		-B {params.output_dir}:/seekdeep_output {input.sif_file} \
		SeekDeep genTargetInfoFromGenomes --primers /input_data/{params.primer_file} \
		--pairedEndLength {params.insert_size} --genomeDir /genome_info/{params.genome_subfolder} \
		--gffDir /genome_info/{params.gff_subfolder} --dout \
		/seekdeep_output/extractedRefSeqs --overWriteDir --numThreads \
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
		sample_names=config['sample_names']
		amino_acid_fnp=config['amino_acid_fnp']
	threads: config['cpus_to_use']
	output:
		setup_done=config['output_folder']+'/finished_setup.txt'
	shell:
		'''
		singularity exec -B {input.data_folder}:/input_data \
		-B {params.output_dir}:/seekdeep_output {input.sif_file} \
		-B {input.genome_root_folder}:/genome_info \
		SeekDeep setupTarAmpAnalysis --samples /input_data/{input.sample_names} \
		--outDir /seekdeep_output/analysis \
		--inputDir /input_data/{params.fastq_folder} \
		--idFile /input_data/{params.primer_file} --lenCutOffs \
		/seekdeep_output/extractedRefSeqs/forSeekDeep/lenCutOffs.txt \
		--overlapStatusFnp /seekdeep_output/extractedRefSeqs/forSeekDeep/refSeqs/overlapStatuses.txt \
		--refSeqsDir /seekdeep_output/extracted_refseqs/forSeekDeep/refSeqs/ \
		--extraExtractorCmds="--checkShortenBars \
		--checkRevComplementForPrimers" \
		--extraQlusterCmds= {config['extra_qluster_cmds']} \
		--extraProcessClusterCmds= {config['extra_process_cluster_cmds']} \
		--numThreads {threads}
		touch {output.setup_done}
		'''

rule runAnalysis:
	input:
		data_folder=config['primer_plus_fastq_binding'],
		sif_file=config['sif_file_location'],
		setup_done=config['output_folder']+'/finished_setup.txt'
	params:
		output_dir=config['output_folder'],
	output:
#		analysis_done=directory(config['output_folder']+'/analysis/popClustering')
		analysis_done=config['output_folder']+'/finished_analysis.txt'
	threads: config['cpus_to_use']
	shell:
		'''
		singularity exec -B {input.data_folder}:/input_data \
		-B {params.output_dir}:/seekdeep_output \
		-H {params.output_dir}/analysis/:/home/analysis \
		{input.sif_file} ./runAnalysis.sh {threads}
		touch {output.analysis_done}
		'''
