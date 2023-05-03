'''
*************************************************************************************************
This is the plasmid_finder_pipe.py script. It uses the PlasmidFinder code from CGE.
It is called by UCT_WGS_pipeline.py.

(c) Dr GJK Marais 2023
**************************************************************************************************
'''

import os

plasmidfinder_path = 'plasmidfinder/plasmidfinder/plasmidfinder.py'
plasmidfinder_db = 'plasmidfinder_db'


def generate_plasmid_finder_commands(assembly_default_dict, output_dir):
    commands = []
    # Iterate through data dictionary and create organism directories
    for org, assemblies in assembly_default_dict.items():
        type_output_directory = os.path.join(output_dir, org)
        if not os.path.exists(type_output_directory):
            os.makedirs(type_output_directory)
        # Iterate through isolates and generates plasmid data for each contig file of the isolate
        for assembly in assemblies:
            org_output_directory = os.path.join(type_output_directory, assembly[0])
            if not os.path.exists(org_output_directory):
                os.makedirs(org_output_directory)
            command = ['python', plasmidfinder_path, '-o', org_output_directory, '-i',
                       assembly[1], '-p',  plasmidfinder_db]
            commands.append((assembly[0], command))
    return commands
