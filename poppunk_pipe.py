'''
*************************************************************************************************
This is the poppunk_pipe.py script. It evaluates the genomic distances between refseq genomes.
It is called by UCT_WGS_pipeline.py.

(c) Dr GJK Marais 2023
**************************************************************************************************
'''

import os
import services


def generate_poppunk_commands(refseq_folders, filename_list, threads):
    commands = []
    for refseq_folder in refseq_folders:
        output_folder_poppunk = os.path.join(refseq_folder[1], 'poppunk_clusters')
        command = ['conda', 'run', '-n', 'poppunk_env', 'poppunk_assign', '--db',
                   services.poppunk_database_folder[refseq_folder[0]], '--query',
                   os.path.join(refseq_folder[1], filename_list), '--output',
                   output_folder_poppunk, '--threads', threads]
        command_distances = ['conda', 'run', '-n', 'poppunk_env', 'poppunk_extract_distances.py', '--distances',
                             os.path.join(output_folder_poppunk, 'poppunk_clusters.dists'), '--output',
                             os.path.join(output_folder_poppunk, 'poppunk_clusters.dists.out')]
        commands.append((command, command_distances))
    return commands
