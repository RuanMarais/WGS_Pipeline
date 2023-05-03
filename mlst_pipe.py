import os
import services

mlst_path = 'mlst_master/mlst.py'
mlst_db = 'mlst_db'


def generate_mlst_commands(assembly_default_dict, output_dir):
    commands = []
    for org, assemblies in assembly_default_dict.items():
        if org in services.org_name_db_mlst:
            type_output_directory = os.path.join(output_dir, org)
            if not os.path.exists(type_output_directory):
                os.makedirs(type_output_directory)
            for assembly in assemblies:
                org_output_directory = os.path.join(type_output_directory, assembly[0])
                if not os.path.exists(org_output_directory):
                    os.makedirs(org_output_directory)
                command = ['python', mlst_path, '-o', org_output_directory, '-i',
                           assembly[1], '-s', services.org_name_db_mlst[org],  '-p', mlst_db, '--tmp_dir',
                           org_output_directory]
                commands.append((assembly[0], command))
    return commands