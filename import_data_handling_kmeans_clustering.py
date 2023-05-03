import csv
from collections import defaultdict


def txt_to_list_of_dicts(file_path_test):
    data = []
    with open(file_path_test, mode='r', encoding='utf-8') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for read_row in reader:
            data.append(read_row)
    return data


def generate_kmeans_data(file_path_poppunk_out):

    # Call the function and store the result in a variable
    list_of_dicts = txt_to_list_of_dicts(file_path_poppunk_out)

    # Convert to array of pairwise distances (core)
    genome_default_dict = defaultdict(list)
    for genome_dict in list_of_dicts:
        core_distance = genome_dict['Core']
        genome = genome_dict['Query']
        genome_default_dict[genome].append(core_distance)

    genome_reference_list = []
    kmeans_data_list = []
    for genome, reference_list in genome_default_dict.items():
        genome_reference_list.append(genome)
        kmeans_data_list.append(reference_list)

    return genome_reference_list, kmeans_data_list, genome_default_dict

