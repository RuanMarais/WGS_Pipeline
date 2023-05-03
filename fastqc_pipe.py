'''
*************************************************************************************************
This is the fastqc_pipe.py script. FastQC is a quality control tool for high throughput sequence data.
A html file is created for each sample.

(c) Dr GJK Marais 2023
**************************************************************************************************
'''
import os


def generate_fastqc_commands(paired_file_list, id_length, output_dir, threads, source_directory):
    commands = []
    # Iterate through trimmed reads
    for paired_files in paired_file_list:
        id_file = paired_files[0][:id_length + 1]
        org_output_directory = os.path.join(output_dir, id_file)
        if not os.path.exists(org_output_directory):
            os.makedirs(org_output_directory)
        fastqc_command = ['fastqc', '-t', str(threads), '-o', org_output_directory,
                          os.path.join(source_directory, paired_files[0]),
                          os.path.join(source_directory, paired_files[1])]
        commands.append((id_file, fastqc_command))
    return commands