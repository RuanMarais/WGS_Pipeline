'''
*************************************************************************************************
This is the quast_pipe.py script. It is used to generate the QUAST quality report for the
assembly of the contigs. It is called by UCT_WGS_pipeline.py.

(c) Dr GJK Marais 2023
**************************************************************************************************
'''

import os
import services


def generate_quast_commands(assembly_default_dict, output_dir, threads, reference_genome_dict):
    commands = []
    # Iterate through data dictionary and create organism directories
    for org, assemblies in assembly_default_dict.items():
        type_output_directory = os.path.join(output_dir, org)
        if not os.path.exists(type_output_directory):
            os.makedirs(type_output_directory)
        # Iterate through isolates and generates quast data for each contig file of the isolate
        for assembly in assemblies:
            org_output_directory = os.path.join(type_output_directory, assembly[0])
            if not os.path.exists(org_output_directory):
                os.makedirs(org_output_directory)
            command = ['conda', 'run', '-n', 'quast_env', 'quast', assembly[1], '-r',
                       reference_genome_dict[org][0], '-o', org_output_directory, '-t', str(threads)]
            commands.append((assembly[0], command))
    return commands