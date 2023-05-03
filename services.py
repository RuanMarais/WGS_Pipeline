import os
from collections import defaultdict
import gzip
import shutil
import re

# org_name_db_resfinder = {'KLEB': "Klebsiella", 'ECLO': "Escherichia coli"}
# reference_genomes_key = {'KLEB': 'reference_genomes/Klebsiella_pneumoniae/CP085866.1.fasta',
#                          'ECLO': 'reference_genomes/Escherichia_coli/ncbi_dataset/data/
#                          GCF_000008865.2/GCF_000008865.2_ASM886v2_genomic.fna',
#                          'CITR': 'reference_genomes/Citrobacter_freundii/ncbi_dataset/data/
#                          GCF_003812345.1/GCF_003812345.1_ASM381234v1_genomic.fna',
#                          'SERR': 'reference_genomes/Serratia_marcescens/ncbi_dataset/data/
#                          GCF_003516165.1/GCF_003516165.1_ASM351616v1_genomic.fna',
#                          'PROV': 'reference_genomes/Providencia_stuartii/ncbi_dataset/data/
#                          GCF_023547145.1/GCF_023547145.1_ASM2354714v1_genomic.fna'}

org_name_key = {'Klebsiella_pneumoniae': 'Klebsiella pneumoniae', 'Escherichia_coli': 'Escherichia coli',
                'Citrobacter_freundii': 'Citrobacter freundii', 'Providencia_stuartii': 'Providencia stuartii',
                'Serratia_marcescens': 'Serratia marcescens'}
org_name_db_mlst = {'Klebsiella_pneumoniae': 'kpneumoniae', 'Escherichia_coli': 'ecoli', 'Citrobacter_freundii': 'cfreundii'}

poppunk_database_folder = {'Klebsiella_pneumoniae': 'poppunk_databases/Klebsiella_pneumoniae_v3_refs', 'Escherichia_coli': 'poppunk_databases/', 'Citrobacter_freundii': 'poppunk_databases/',
                           'Providencia_stuartii': 'poppunk_databases/', 'Serratia_marcescens': 'poppunk_databases/', 'Acinetobacter_baumanii': 'poppunk_databases/',
                           'Mycobacterium_tuberculosis': 'poppunk_databases/'}

quast_keys = ['Assembly', '# contigs (>= 0 bp)', '# contigs (>= 1000 bp)', '# contigs (>= 5000 bp)',
              '# contigs (>= 10000 bp)', '# contigs (>= 25000 bp)', '# contigs (>= 50000 bp)',
              'Total length (>= 0 bp)', 'Total length (>= 1000 bp)', 'Total length (>= 5000 bp)',
              'Total length (>= 10000 bp)', 'Total length (>= 25000 bp)', 'Total length (>= 50000 bp)',
              '# contigs', 'Largest contig', 'Total length', 'GC (%)', 'N50', 'NG50', 'N75', 'NG75', 'L50', 'LG50',
              'L75', 'LG75', '# misassemblies', '# misassemblies contigs', 'Misassembled contigs length',
              '# local misassemblies', '# scaffold gap ext. mis.', '# scaffold gap loc. mis.',
              '# unaligned mis. contigs', '# unaligned contigs', 'Unaligned length', 'Genome fraction (%)',
              'Duplication ratio', '# N\'s per 100 kbp', '# mismatches per 100 kbp', '# indels per 100 kbp',
              'Largest alignment', 'Total aligned length', 'NA50', 'NGA50', 'NA75', 'NGA75',
              'LA50', 'LGA50', 'LA75', 'LGA75']

vfdb_available = ['Klebsiella_pneumoniae', 'Escherichia_coli']


def filename_list_generate(file_identifier, folder_path):
    filenames_output = []
    for file in os.listdir(folder_path):
        if os.path.isfile(os.path.join(folder_path, file)):
            if file_identifier in file:
                filenames_output.append(file)
    return filenames_output


def paired_read_list_generate(id_length, file_id_r1, file_id_r2, filename_list):
    # Separate filenames by id
    sample_sorted_dict = defaultdict(list)
    for file in filename_list:
        unique_filename = file[:id_length + 1]
        sample_sorted_dict[unique_filename].append(file)

    # Separate pairs
    paired_sorted_list = []
    for key, files in sample_sorted_dict.items():
        read_1 = None
        read_2 = None
        for file in files:
            if file_id_r1 in file:
                read_1 = file
            elif file_id_r2 in file:
                read_2 = file
        if read_1 is not None and read_2 is not None:
            output = (read_1, read_2)
            paired_sorted_list.append(output)

    return paired_sorted_list


def generate_result_info_dict(result_folder_path):
    result_data_structure = defaultdict(list)
    folders = os.listdir(result_folder_path)
    directory_paths = [(os.path.join(result_folder_path, folder), folder)
                       for folder in folders
                       if os.path.isdir(os.path.join(result_folder_path, folder))]
    for directory in directory_paths:
        results = os.listdir(directory[0])
        result_data_structure[directory[1]].extend([(os.path.join(directory[0], result), result)
                                                    for result in results
                                                    if os.path.isdir(os.path.join(directory[0], result))])
    return result_data_structure


def decompress_gzip_file(gzip_path, dest_path, error_log):
    try:
        with gzip.open(gzip_path, 'rb') as gzip_file:
            with open(dest_path, 'wb') as decompressed_file:
                shutil.copyfileobj(gzip_file, decompressed_file)
    except:
        error_log.append(f'Error decompressing file: {gzip_path}')


# This formats filename to replace all special characters with underscores
def format_filename_to_dashes(filename):
    # Replace dashes, periods, and other non-alphanumeric characters with underscores
    return re.sub(r'[^a-zA-Z0-9]', '_', filename)
