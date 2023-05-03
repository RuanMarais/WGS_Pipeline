'''
*************************************************************************************************
This is the spades_pipe.py script. This script generates the spades commands for each sample.

(c) Dr GJK Marais 2023
**************************************************************************************************
'''

import os


def generate_spades_commands(paired_file_list, id_length, threads, output_directory, source_directory):
    output = []
    for paired_files in paired_file_list:
        id_file = paired_files[0][:id_length + 1]
        result_file = os.path.join(output_directory, id_file)
        spades_command = ['spades.py', '-t', str(threads), '--careful', '-1', os.path.join(source_directory, paired_files[0]), '-2',
                          os.path.join(source_directory, paired_files[1]), '-o', result_file]
        output.append((id_file, spades_command))
    return output
