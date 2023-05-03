'''
*************************************************************************************************
This is the roary_plot_pipe.py script. It processes Roary and Fasttree generated files for figure
output.
It is called by UCT_WGS_pipeline.py.

(c) Dr GJK Marais 2023
**************************************************************************************************
'''


def generate_roary_plot_commands(roary_dict, fasttree_output_dir):
    commands = []
    for org, analysis in roary_dict.items():
        command_one = ['python', 'roary_plots.py', '--labels', '--format', 'svg', f'{fasttree_output_dir}/{org}.newick',
                       analysis[0]]
        command_two = ['python', 'roary_plots.py', '--labels', f'{fasttree_output_dir}/{org}.newick', analysis[0]]
        commands.append((command_one, command_two))
    return commands
