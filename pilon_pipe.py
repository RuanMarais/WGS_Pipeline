'''
*************************************************************************************************
This is the pilon_pipe.py script. It is used to polish the assembly.
It is called by UCT_WGS_pipeline.py.

(c) Dr GJK Marais 2023
**************************************************************************************************
'''

import os


def align_bowtie(assembly, forward_read, reverse_read, org_output_directory, threads, org_id):
    sam_file = os.path.join(org_output_directory, org_id + '.sam')
    bam_file = os.path.join(org_output_directory, org_id + '.bam')
    command_build = ['bowtie2-build', assembly, assembly]
    command_align = ['bowtie2', '--local', '-p', str(threads), '-x', assembly, '-1', forward_read, '-2', reverse_read, '-S', sam_file]
    command_samtools_sort = ['samtools', 'sort', sam_file, '-o', bam_file]
    command_samtools_index = ['samtools', 'index', bam_file]
    command_samtools_rm = ['rm', sam_file]
    return [command_build, command_align, command_samtools_sort, command_samtools_index, command_samtools_rm]


def generate_pilon_commands(assemblies, output_dir, threads, general_outdir):
    commands = []
    # Iterate through data dictionary and create organism directories
    for assembly in assemblies:
        # Iterate through isolates and generates pilon commands for each isolate assembly file
        org_output_directory = os.path.join(output_dir, assembly[0])
        if not os.path.exists(org_output_directory):
            os.makedirs(org_output_directory)
        align_commands = align_bowtie(assembly[1],  f'{general_outdir}/trim_results/' + assembly[0] + '_paired_1.fastq',
                                      f'{general_outdir}/trim_results/' + assembly[0] + '_paired_2.fastq',
                                      org_output_directory, threads, assembly[0])
        command_pilon = ['pilon', '--genome', assembly[1], '--bam',
                         os.path.join(org_output_directory, assembly[0] + '.bam'), '--changes',
                         '--output', assembly[0], '--outdir', org_output_directory, '--threads', str(threads)]
        commands_out = [align_commands, command_pilon]
        commands.append([assembly[0], commands_out])
    return commands
