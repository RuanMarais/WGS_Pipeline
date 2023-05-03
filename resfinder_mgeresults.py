'''
*************************************************************************************************
This is the resfinder_mgeresults.py script. Its uses the ResFinder tool from CGE to detect
resistance genes in mobile genetic elements.
It is called by UCT_WGS_pipeline.py.

(c) Dr GJK Marais 2023
**************************************************************************************************
'''

import os

resfinder_path = 'resfinder_master/src/resfinder/run_resfinder.py'
resfinder_db_path = 'resfinder_db'


def generate_resfinder_commands(mge_results_dict, output_dir):
    commands = []
    for org, sequence_file in mge_results_dict.items():
        type_output_directory = os.path.join(output_dir, org)
        if not os.path.exists(type_output_directory):
            os.makedirs(type_output_directory)
        command = ['python', resfinder_path, '-o', type_output_directory, '-l', str(0.6), '-t', str(0.9),
                   '--acquired', '-ifa', sequence_file[0], '-db_res', resfinder_db_path]
        commands.append(command)
    return commands