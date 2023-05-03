import Bio.SeqIO as SeqIO
from collections import defaultdict
import os
import pandas as pd
import subprocess
import shutil

# test_paths = ['parsnp_data/klep_2018to2022_complete_good.fasta']
# test_output = 'parsnp_data/kpneumoniae_tree_data'

# def parse_bvbrc_sequences_into_separate_fasta_files(filepaths, output_dir):
#     for filepath in filepaths:
#         with open(filepath, 'r') as f:
#             sequences = SeqIO.parse(f, 'fasta')
#             unique_sequences = defaultdict(list)
#             for sequence in sequences:
#                 header = sequence.description
#                 bracket_index = header.find('[')
#                 bracket_index_end = header.find(']')
#                 seq_id = header[bracket_index+1:bracket_index_end]
#                 unique_sequences[seq_id].append(sequence)
#
#             for key, sequences in unique_sequences.items():
#                 with open(os.path.join(output_dir, f'{key}.fasta'), 'x') as f:
#                     SeqIO.write(sequences, f, 'fasta')

def generate_evaluation_folder_parsnp(assembly_default_dict, output_dir):
    for org, assemblies in assembly_default_dict.items():
        type_output_directory = os.path.join(output_dir, org)
        if not os.path.exists(type_output_directory):
            os.makedirs(type_output_directory)
        for assembly in assemblies:
            shutil.copy(assembly[1], os.path.join(type_output_directory, assembly[0] + '.fasta'))
