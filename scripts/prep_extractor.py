import subprocess
extractor_commands=snakemake.input['extractor_commands']
output_folder=snakemake.params['output_folder']
subprocess.call(f'mkdir -p {output_folder}', shell=True)
for line in open(extractor_commands):
	print(line)
	split_line=line.strip().split()
	sample_index=split_line.index('--sampleName')+1
	sample=split_line[sample_index][3:]
	output_file=open(output_folder+f'/{sample}_extraction_command.sh', 'w')
	output_file.write(line)
	print(sample)
