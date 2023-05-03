'''
*************************************************************************************************
This is the kaptive_pipe.py script. Its uses the Kaptive tool to detect capsule type.
It is called by UCT_WGS_pipeline.py.

(c) Dr GJK Marais 2023
**************************************************************************************************
'''

import os
import std_paths


def generate_kaptive_commands(assembly_default_dict, output_dir):
    commands = []
    for org, assemblies in assembly_default_dict.items():
        if org in std_paths.kaptive_db_path:
            type_output_directory = os.path.join(output_dir, org)
            if not os.path.exists(type_output_directory):
                os.makedirs(type_output_directory)
            for assembly in assemblies:
                org_output_directory = os.path.join(type_output_directory, assembly[0])
                if not os.path.exists(org_output_directory):
                    os.makedirs(org_output_directory)
                for kapsular_type_db in std_paths.kaptive_db_path[org]:
                    prefix = kapsular_type_db[0]
                    command = ['python', std_paths.kaptive_path, '-a', assembly[1], '-k',
                               kapsular_type_db[1], '--threads', '1', '-o',
                               os.path.join(org_output_directory, prefix)]
                    commands.append((assembly[0], command))
    return commands