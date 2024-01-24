import subprocess
process_cluster_commands=snakemake.input['process_cluster_commands']
output_folder=snakemake.params['output_folder']
subprocess.call(f'mkdir -p {output_folder}', shell=True)
analysis_dir=snakemake.params.analysis_dir
for line in open(process_cluster_commands):
	print(line)
	split_line=line.strip().split()
	sample_index=split_line.index('--experimentName')+1
	sample=split_line[sample_index]
	output_file=open(output_folder+f'/{sample}_process_cluster_command.sh', 'w')
	output_file.write('cd {analysis_dir}\n')
	output_file.write(line)
	print(sample)
