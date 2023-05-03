'''
*************************************************************************************************
This is the kmerfinder_pipe.py script. It uses the KmerFinder code from CGE to identify the species.
It is called by UCT_WGS_pipeline.py.

(c) Dr GJK Marais 2023
**************************************************************************************************
'''

import os
import std_paths


def generate_kmerfinder_commands(assembly_default_dict, output_dir):
    commands = []
    for org, assemblies in assembly_default_dict.items():
        type_output_directory = os.path.join(output_dir, org)
        if not os.path.exists(type_output_directory):
            os.makedirs(type_output_directory)
        for assembly in assemblies:
            org_output_directory = os.path.join(type_output_directory, assembly[0])
            if not os.path.exists(org_output_directory):
                os.makedirs(org_output_directory)
            command = ['python', std_paths.kmerfinder_path, '-i', assembly[1], '-o', org_output_directory,
                       '-db', std_paths.kmerfinder_db_path, '-tax', std_paths.kmerfinder_db_tax_path, '-x']
            commands.append((assembly[0], command))
    return commands