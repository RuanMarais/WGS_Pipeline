import os
import subprocess


def generate_bwa_index(candidate_reference_defaultdict):
    commands = []
    for org, data in candidate_reference_defaultdict.items():
        for item in data:
            bwa_index_command = ['bwa', 'index', item[1]]
            commands.append(bwa_index_command)
    return commands


def align_scaffolds_bwa(data_dict, threads, output_directory):
    commands = []
    for org, data in data_dict.items():
        references = data[0]
        assemblies = data[1]
        type_output_directory = os.path.join(output_directory, org)
        if not os.path.exists(type_output_directory):
            os.makedirs(type_output_directory)
        for reference in references:
            reference_output_directory = os.path.join(type_output_directory, reference[2])
            if not os.path.exists(reference_output_directory):
                os.makedirs(reference_output_directory)
            for assembly in assemblies:
                output_sam = os.path.join(reference_output_directory, f'{assembly[0]}.sam')
                bwa_command = ['bwa', 'mem', '-t', str(threads), reference[1], assembly[1], '-o', output_sam]
                convert_to_bam = ['samtools', 'view', '-b', '-S', output_sam]
                output_bam_file = os.path.join(reference_output_directory, f'{assembly[0]}.bam')
                sorted_output_bam_file = os.path.join(reference_output_directory, f'{assembly[0]}.sorted.bam')
                sort_bam = ['samtools', 'sort', output_bam_file, '-o', sorted_output_bam_file]
                index_bam = ['samtools', 'index', sorted_output_bam_file]
                # Call variants
                mpile_up = ['bcftools', 'mpileup', '-Ou', '-f', reference[1], sorted_output_bam_file]
                vcf_file = os.path.join(reference_output_directory, f'{assembly[0]}.vcf.gz')
                call_process = ['bcftools', 'call', '--ploidy', '1', '-mv', '-Oz', '-o', vcf_file]
                index_call = ['bcftools', 'index', vcf_file]
                # Normalise indels
                bcf_file = os.path.join(reference_output_directory, f'{assembly[0]}.norm.bcf')
                bcf_file_filtered = os.path.join(reference_output_directory, f'{assembly[0]}.norm.flt-indels.bcf')
                norm_indel = ['bcftools', 'norm', '-f', reference[1], vcf_file, '-Ob', '-o', bcf_file]
                bcf_filter = ['bcftools', 'filter', '--IndelGap', '5', bcf_file, '-Ob', '-o', bcf_file_filtered]
                # Generate consensus sequence
                cat_ref = ['cat', reference[1]]
                generate_consensus = ['bcftools', 'consensus', vcf_file]
                consensus_file = os.path.join(reference_output_directory, f'{assembly[0]}.consensus.fasta')
                mummer = ['dnadiff', reference[1], consensus_file, '-p', os.path.join(reference_output_directory,
                                                                                      f'{assembly[0]}.mummer')]
                commands_to_add = {'BWA': bwa_command, 'toBAM': convert_to_bam, 'BAM': output_bam_file,
                                   'sortBAM': sort_bam, 'indexBAM': index_bam, 'mpilup': mpile_up,
                                   'indexVCF': index_call, 'normIndel': norm_indel, 'bcfFilter': bcf_filter,
                                   'VCF': vcf_file, 'genVCF': call_process, 'consensusfile': consensus_file,
                                   'genConsensus': generate_consensus, 'catRef': cat_ref, 'mummer': mummer}
                commands.append(commands_to_add)
    return commands


def generate_consensus_sequences(subprocess_commands_dict_list, error_log):
    for command_dict in subprocess_commands_dict_list:
        try:
            subprocess.run(command_dict['BWA'], check=True)
        except subprocess.CalledProcessError as e:
            error_log.append(f"Candidate reference index for BWA: {command_dict['BWA']} with error code: {e.returncode}")
        try:
            subprocess.run(command_dict['toBAM'], check=True, stdout=open(command_dict['BAM'], 'wb'))
        except subprocess.CalledProcessError as e:
            error_log.append(f"Conversion to BAM: {command_dict['toBAM']} with error code: {e.returncode}")
        try:
            subprocess.run(command_dict['sortBAM'], check=True)
        except subprocess.CalledProcessError as e:
            error_log.append(f"Sorting BAM: {command_dict['sortBAM']} with error code: {e.returncode}")
        try:
            subprocess.run(command_dict['indexBAM'], check=True)
        except subprocess.CalledProcessError as e:
            error_log.append(f"Indexing BAM: {command_dict['indexBAM']} with error code: {e.returncode}")
        try:
            with subprocess.Popen(command_dict['mpilup'], stdout=subprocess.PIPE) as proc1:
                subprocess.run(command_dict['genVCF'], stdin=proc1.stdout)
        except:
            error_log.append('Generating VCF failed')
        try:
            subprocess.run(command_dict['indexVCF'], check=True)
        except subprocess.CalledProcessError as e:
            error_log.append(f"Indexing VCF: {command_dict['indexVCF']} with error code: {e.returncode}")
        try:
            subprocess.run(command_dict['normIndel'], check=True)
        except subprocess.CalledProcessError as e:
            error_log.append(f"Normalisation of Indels failed: {command_dict['normIndel']} with error code: {e.returncode}")
        try:
            subprocess.run(command_dict['bcfFilter'], check=True)
        except subprocess.CalledProcessError as e:
            error_log.append(f"Filtering BCF failed: {command_dict['bcfFilter']} with error code: {e.returncode}")
        try:
            subprocess.run(command_dict['indexVCF'], check=True)
        except subprocess.CalledProcessError as e:
            error_log.append(f"Indexing VCF: {command_dict['indexVCF']} with error code: {e.returncode}")
        try:
            with subprocess.Popen(command_dict['catRef'], stdout=subprocess.PIPE) as proc1:
                with subprocess.Popen(command_dict['genConsensus'], stdin=proc1.stdout, stdout=subprocess.PIPE) as proc2:
                    with open(command_dict['consensusfile'], 'wb') as f_out:
                        f_out.write(proc2.stdout.read())
        except:
            error_log.append('Consensus generation failed')
        try:
            subprocess.run(command_dict['mummer'], check=True)
        except subprocess.CalledProcessError as e:
            error_log.append(f"Mummer command: {command_dict['mummer']} failed with error code: {e.returncode}")


def find_optimal_reference():
    pass