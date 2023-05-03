'''
*************************************************************************************************
This is the resfinder_pipe.py script. Its uses the ResFinder tool from CGE to detect resistance genes.
It is called by UCT_WGS_pipeline.py.

(c) Dr GJK Marais 2023
**************************************************************************************************
'''

import os

resfinder_path = 'resfinder_master/src/resfinder/run_resfinder.py'
resfinder_db_path = 'resfinder_db'
disinfinder_db_path = 'disinfinder_db'
pointfinder_db_path = 'pointfinder_db'


def generate_resfinder_commands(assembly_default_dict, output_dir):
    commands = []
    for org, assemblies in assembly_default_dict.items():
        type_output_directory = os.path.join(output_dir, org)
        if not os.path.exists(type_output_directory):
            os.makedirs(type_output_directory)
        for assembly in assemblies:
            org_output_directory = os.path.join(type_output_directory, assembly[0])
            if not os.path.exists(org_output_directory):
                os.makedirs(org_output_directory)
            command = ['python', resfinder_path, '-o', org_output_directory, '-l', str(0.6), '-t', str(0.9),
                       '--acquired', '-ifa', assembly[1], '-db_res', resfinder_db_path, '-db_disinf',
                       disinfinder_db_path, '-db_point', pointfinder_db_path]
            commands.append((assembly[0], command))
    return commands
