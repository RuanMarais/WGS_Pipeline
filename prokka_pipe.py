'''
*************************************************************************************************
This is the prokka_pipe.py script. It uses the Prokka code to annotate resistance gene containing
contigs. It is called by UCT_WGS_pipeline.py.

(c) Dr GJK Marais 2023
**************************************************************************************************
'''
import os


def generate_prokka_commands(assembly_default_dict, output_dir):
    commands = []
    # Iterate through data dictionary and create organism directories
    for org, assemblies in assembly_default_dict.items():
        type_output_directory = os.path.join(output_dir, org)
        if not os.path.exists(type_output_directory):
            os.makedirs(type_output_directory)
        # Iterate through isolates and generates prokka annotation files for each contig file
        for assembly in assemblies:
            out_path = os.path.join(type_output_directory, assembly[0])
            genus_species = org.split('_')
            command = ['conda', 'run', '-n', 'prokka_env', 'prokka', '--outdir', out_path, '--compliant',
                       '--centre', 'UCT_Micro', '--genus', genus_species[0], '--species', genus_species[1],
                       '--prefix', assembly[0], assembly[1]]
            commands.append((assembly[0], command))
    return commands

