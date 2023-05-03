'''
*************************************************************************************************
This is the hts_stream.py script. This script generates the hts stream commands including
de-duplication.

(c) Dr GJK Marais 2023
**************************************************************************************************
'''
import os

def generate_deduplicator_commands(paired_file_list, id_length, folder_path, output_dir):
    commands = []
    # Iterate through trimmed reads
    for paired_files in paired_file_list:
        id_file = paired_files[0][:id_length + 1]
        dedup_command = ['hts_SuperDeduper', '-1', os.path.join(folder_path, paired_files[0]), '-2',
                         os.path.join(folder_path, paired_files[1]),  '-f', os.path.join(output_dir, id_file),
                         '-A', os.path.join(output_dir, 'deduplicator_log')]
        commands.append(dedup_command)
    return commands