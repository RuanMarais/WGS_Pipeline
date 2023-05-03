'''
*************************************************************************************************
This is the gen_find.py script. It uses the CGE mydbfinder code to query provided databases

(c) Dr GJK Marais 2023
**************************************************************************************************
'''
import os


gen_find_path = 'general_dbfinder_master/mydbfinder/mydbfinder.py'


def generate_gen_db_commands(assembly_default_dict, db_path, db, output_dir):
    commands = []
    # Iterate through data dictionary and create organism directories
    for org, assemblies in assembly_default_dict.items():
        type_output_directory = os.path.join(output_dir, org)
        if not os.path.exists(type_output_directory):
            os.makedirs(type_output_directory)
        # Iterate through isolates and generates evaluates for virulence factors
        for assembly in assemblies:
            org_output_directory = os.path.join(type_output_directory, assembly[0])
            if not os.path.exists(org_output_directory):
                os.makedirs(org_output_directory)
            command = ['python', gen_find_path, '-i', assembly[1], '-o', org_output_directory, '-p', db_path,
                       '-d', db, '-l', '0.9', '-t', '0.98']
            commands.append((assembly[0], command))
    return commands