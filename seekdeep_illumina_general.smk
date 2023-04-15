configfile: 'seekdeep_nanopore_general.yaml'

rule all:
	input:
#		primer_info=config['output_folder']+'/extractedRefSeqs/locationsByGenome/Pf3D7_infos.tab.txt'
#		done=config['output_folder']+'/finished_setup.txt'
		analysis_done=config['output_folder']+'/finished_analysis.txt'	
#		analysis_done=config['output_folder']+'/analysis/popClustering'

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
		sif_file=config['sif_file_location'],
		primer_info=config['output_folder']+'/extractedRefSeqs/locationsByGenome/Pf3D7_infos.tab.txt'
	params:
		output_dir=config['output_folder'],
		primer_file=config['primer_file'],
		fastq_folder=config['fastq_subfolder']
	output:
		setup_done=config['output_folder']+'/finished_setup.txt'
	params:
		output_dir=config['output_folder'],
		primer_file=config['primer_file'],
		fastq_folder=config['fastq_subfolder'],
		fake_insert_size=config['fake_insert_size']
#		done='output_files/setupTarAmpAnalysis_finished.txt'
	threads: config['cpus_to_use']
	shell:
		'''
		singularity exec -B {input.data_folder}:/input_data \
		-B {params.output_dir}:/seekdeep_output {input.sif_file} \
		SeekDeep setupTarAmpAnalysis  --outDir \
		/seekdeep_output/analysis --technology nanopore \
		--uniqueKmersPerTarget \
		/seekdeep_output/extractedRefSeqs/forSeekDeep/uniqueKmers.tab.txt.gz \
		--idFile /input_data/{params.primer_file} --lenCutOffs \
		/seekdeep_output/extractedRefSeqs/forSeekDeep/lenCutOffs.txt \
		--inputDir /input_data/{params.fastq_folder} \
		--doNotGuessRecFlags --extraExtractorCmds="--primerWithinStart 40 \
		--minLenCutOff 100" --extraKlusterCmds="--readLengthMinDiff 30 \
		--numThreads {threads}" \
		--extraProcessClusterCmds="--replicateMinTotalReadCutOff 10 \
		--previousPop-largeBaseIndel 0.99 --previousPop-oneBaseIndel 0.99 \
		--previousPop-twoBaseIndel 0.99" \
		--previousPopSeqsDir /seekdeep_output/extractedRefSeqs/forSeekDeep/refSeqs/ \
		--numThreads {threads}
		touch {output.setup_done}
		'''
todo: merge above with below:
		'''
		SeekDeep setupTarAmpAnalysis --samples {input.sample_names} \
		--outDir {output.analysis_dir} \
		--inputDir {input.fastq_dir} \
		--idFile {input.id_file} \
		--lenCutOffs {input.len_cutoffs} \
		--overlapStatusFnp {input.overlap_stats} \
		--refSeqsDir {input.refseqs} \
		--extraExtractorCmds="--checkShortenBars \
		--checkRevComplementForPrimers" \
		--extraQlusterCmds="--illuminaAllowHomopolyers" \
		--extraProcessClusterCmds="--lowFreqHaplotypeFracCutOff 0.01 \
		--gffFnp {input.genome_folder}/info/gff/Pf3D7.gff \
		--genomeFnp {input.genome_folder}/genomes/Pf3D7.fasta \
		--knownAminoAcidChangesFnp {input.genome_folder}/info/pf_drug_resistant_aaPositions_k13_updated.tsv"
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
		{input.sif_file} ./runAnalysis.sh
		touch {output.analysis_done}
		'''

#		echo "cd /seekdeep_output/analysis\n./runAnalysis.sh" > {params.output_dir}/analysis_runner.sh
#		chmod +x {params.output_dir}/analysis_runner.sh
#		singularity exec -B {input.data_folder}:/input_data \
#		-B {params.output_dir}:/seekdeep_output {input.sif_file} \
#		/seekdeep_output/analysis_runner.sh {threads}
#		touch {output.analysis_done}

