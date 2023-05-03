'''
*************************************************************************************************
This is the ragtag_pipe.py script. Its uses the Ragtag tool to perform reference guided scaffolding.
It is called by UCT_WGS_pipeline.py.

(c) Dr GJK Marais 2023
**************************************************************************************************
'''

import os

import services


def generate_ragtag_commands(assemblies_dict, output_dir, threads, reference_dict):
    commands = []
    for org, assemblies in assemblies_dict.items():
        for assembly in assemblies:
            org_output_directory = os.path.join(output_dir, assembly[0])
            if not os.path.exists(org_output_directory):
                os.makedirs(org_output_directory)
            command = ['conda', 'run', '-n', 'ragtag_env', 'ragtag.py', 'scaffold', reference_dict[org][0], assembly[1],
                       '-o', org_output_directory, '-t', str(threads)]
            commands.append((assembly[0], command))
    return commands