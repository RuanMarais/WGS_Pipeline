'''
*************************************************************************************************
This is the roary_pipe.py script. Its uses the roary tool to perform pangenome analysis on prokka
produced results.
It is called by UCT_WGS_pipeline.py.

(c) Dr GJK Marais 2023
**************************************************************************************************
'''
import os
import glob

def generate_roary_commands(analysis_folder, threads):
    commands = []
    folders = os.listdir(analysis_folder)
    directory_paths = [(os.path.join(analysis_folder, folder), folder) for folder in folders if
                       os.path.isdir(os.path.join(analysis_folder, folder))]
    for directory in directory_paths:
        files_pass = glob.glob(f'{directory[0]}/*.gff')
        command = ['roary', '-p', str(threads), '-f', directory[0], '-e', '-n', '--mafft', *files_pass]
        commands.append(command)
    return commands