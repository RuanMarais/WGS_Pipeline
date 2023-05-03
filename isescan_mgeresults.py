'''
*************************************************************************************************
This is the isescan_mgeresults.py script. This script annotates insertion sequences.

(c) Dr GJK Marais 2023
**************************************************************************************************
'''


import os


def generate_isescan_commands(is_sequences_dict, output_dir, threads):
    commands = []
    # Iterate through IS sequences
    for org, sequence_file in is_sequences_dict.items():
        type_directory = os.path.join(output_dir, org)
        if not os.path.exists(type_directory):
            os.makedirs(type_directory)
        command = ['conda', 'run', '-n', 'isescan_env', 'isescan.py', '--seqfile', sequence_file[0],
                   '--output', type_directory, '--nthread', str(threads)]
        commands.append(command)
    return commands
