'''
*************************************************************************************************
This is the mgefinder_pipe.py script. This script generates the commands required for mgefinder.

(c) Dr GJK Marais 2023
**************************************************************************************************
'''

import os
import shutil
import services


def generate_bwa_commands(paired_file_list, id_length, source_folder, output_dir, threads, org_name_dict,
                          reference_location_dict):
    commands = []
    index_list = []
    mge_commands = []
    reference_location_dict_mge = {}
    # Iterate through trimmed reads
    for paired_files in paired_file_list:
        id_file = paired_files[0][:id_length + 1]
        org_name = org_name_dict[id_file]

        # Create mgefinder source directories
        type_output_directory = os.path.join(output_dir, org_name)
        if not os.path.exists(type_output_directory):
            os.makedirs(type_output_directory)
        processing_folder_assembly = os.path.join(type_output_directory, '00.assembly')
        if not os.path.exists(processing_folder_assembly):
            os.makedirs(processing_folder_assembly)
        processing_folder_bam = os.path.join(type_output_directory, '00.bam')
        if not os.path.exists(processing_folder_bam):
            os.makedirs(processing_folder_bam)
        processing_folder_genome = os.path.join(type_output_directory, '00.genome')
        if not os.path.exists(processing_folder_genome):
            os.makedirs(processing_folder_genome)

        reference_location_mge = os.path.join(processing_folder_genome, f'{org_name}.fasta')
        reference_location_dict_mge[org_name] = reference_location_mge
        if not os.listdir(processing_folder_genome):
            shutil.copy(reference_location_dict[org_name][0], reference_location_mge)

        bwa_index_command = ['bwa', 'index', reference_location_dict_mge[org_name]]
        bwa_command = ['bwa', 'mem', '-t', str(threads), reference_location_dict_mge[org_name],
                       os.path.join(source_folder, paired_files[0]), os.path.join(source_folder, paired_files[1]), '-o',
                       os.path.join(processing_folder_bam, f'{id_file}.{org_name}.sam')]
        format_command = ['conda', 'run', '-n', 'mgefinder', 'mgefinder', 'formatbam',
                          os.path.join(processing_folder_bam, f'{id_file}.{org_name}.sam'),
                          os.path.join(processing_folder_bam, f'{id_file}.{org_name}.bam')]
        mge_run_command = ['conda', 'run', '-n', 'mgefinder', 'mgefinder', 'workflow', 'denovo',
                           '--cores', str(threads), type_output_directory]
        if bwa_index_command not in index_list:
            index_list.append(bwa_index_command)
        if mge_run_command not in mge_commands:
            mge_commands.append(mge_run_command)
        commands.append([bwa_command, format_command])
    return index_list, commands, mge_commands


# def analyse_mgefinder_output(output_dir):
#     mgefinder_output = []
#     for org_name in services.reference_genomes_key:
#         type_output_directory = os.path.join(output_dir, org_name)
#         processing_folder_bam = os.path.join(type_output_directory, '00.bam')
#         for file in os.listdir(processing_folder_bam):
#             if file.endswith('.bam'):
#                 id_file = file.split('.')[0]
#                 mgefinder_output.append([id_file, org_name, os.path.join(processing_folder_bam, file)])
#     return mgefinder_output