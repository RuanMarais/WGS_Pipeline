'''
*************************************************************************************************
This is the ncbi_refseq_download_pipe.py script. It retrieves all complete refseq genomes from NCBI.
It is called by UCT_WGS_pipeline.py.

(c) Dr GJK Marais 2023
**************************************************************************************************
'''

import os
import shlex


def download_genomes(output_dir, orgs_to_download):
    commands = []
    for org in orgs_to_download:
        type_output_directory = os.path.join(output_dir, org)
        input_org = org.replace('_', ' ')
        if not os.path.exists(type_output_directory):
            os.makedirs(type_output_directory)
        command = f'ncbi-genome-download --formats fasta --assembly-levels complete -o {type_output_directory} -v --flat-output --genera "{input_org}" bacteria'
        commands.append(command)
    return commands