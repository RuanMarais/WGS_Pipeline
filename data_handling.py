import csv
from collections import defaultdict
import json
import os
import pandas as pd
import services
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import virulence_data_handling
import std_paths
import networkx as nx


def generate_quast_datadicts(data_file):
    """
    Generate the data dictionaries from the raw output of QUAST with keys either representing
    the analysis type
    :param data_file:
    :return:
    """
    output_dict = {}
    try:
        with open(data_file, 'r') as tsvfile:
            reader = csv.DictReader(tsvfile, delimiter='\t')
            output_dict = [res_dict for res_dict in reader]
    except:
        with open('test/error_log', 'a+') as f:
            f.seek(0, 2)
            f.write(f'Generate quast datadicts failed for: {data_file}\n')
    return output_dict[0]


def generate_kaptive_datadicts(data_file):
    output_kaptive = None
    try:
        with open(data_file, 'r') as f:
            json_string = f.read()
            data = json.loads(json_string)
            output_kaptive = (data[0]['Best match']['Type'], data[0]['Best match']['Match confidence'])
    except:
        with open('test/error_log', 'a+') as f:
            f.seek(0, 2)
            f.write(f'Generate kaptive datadicts failed for: {data_file}\n')
    return output_kaptive


def generate_resistance_datadicts(data_file):
    """
    Generate the data dictionaries from the raw output of ResFinder with keys either
    representing resistance genes or contigs
    :param data_file:
    :return:
    """
    output_dict_gene_sorted = defaultdict(list)
    output_dict_contig_sorted = defaultdict(list)
    try:
        with open(data_file, 'r') as txtfile:
            reader = csv.DictReader(txtfile, delimiter='\t')
            raw_data = [res_dict for res_dict in reader]
            for item in raw_data:
                if 'Resistance gene' in item:
                    output_dict_gene_sorted[item['Resistance gene']].append(item)
                if 'Contig' in item:
                    output_dict_contig_sorted[item['Contig']].append(item['Resistance gene'])
    except:
        with open('test/error_log', 'a+') as f:
            f.seek(0, 2)
            f.write(f'Generate resistance datadicts failed for: {data_file}\n')
    return output_dict_gene_sorted, output_dict_contig_sorted


def generate_plasmid_datadicts(data_file):

    output_dict_plasmid_sorted = defaultdict(list)
    output_dict_contig_sorted = defaultdict(list)
    try:
        with open(data_file, 'r') as f:
            json_string = f.read()
            data = json.loads(json_string)
        for key, vals in data['plasmidfinder']['results']['Enterobacteriales']['enterobacteriales'].items():
            output_dict_plasmid_sorted[vals['plasmid']].append(vals['contig_name'])
            output_dict_contig_sorted[vals['contig_name']].append(vals['plasmid'])
    except:
        with open('test/error_log', 'a+') as f:
            f.seek(0, 2)
            f.write(f'Generate plasmid datadicts failed for: {data_file}\n')
    return output_dict_plasmid_sorted, output_dict_contig_sorted


def generate_mlst_datadicts(data_file):
    output_mlst = None
    try:
        with open(data_file, 'r') as f:
            json_string = f.read()
            data = json.loads(json_string)
            output_mlst = (data['mlst']['user_input']['organism'], data['mlst']['results']['sequence_type'])
    except:
        with open('test/error_log', 'a+') as f:
            f.seek(0, 2)
            f.write(f'Generate mlst datadicts failed for: {data_file}\n')
    return output_mlst


def generate_heat_map(image_name, x_axis_name, default_dict_data_item, output_directory):
    df_gene_list = []
    isolates = set()
    heat_map_data_dictionary = {}
    for gene, isolate_list in default_dict_data_item.items():
        df_gene_list.append(gene)
        for isolate in isolate_list:
            isolates.add(isolate)
    df_gene_list_sorted = sorted(df_gene_list)
    for isolate in isolates:
        gene_profile = []
        for gene in df_gene_list_sorted:
            if isolate in default_dict_data_item[gene]:
                gene_profile.append(1)
            else:
                gene_profile.append(0)
        heat_map_data_dictionary[isolate] = gene_profile
    
    if len(df_gene_list) > 1 and len(isolates) > 1:
        sns.set(font_scale=0.3)
        df = pd.DataFrame(heat_map_data_dictionary, index=df_gene_list_sorted).transpose()
        ax = sns.clustermap(df, linewidth=.5, linecolor='gray', col_cluster=False, row_cluster=True,
                            metric='jaccard', xticklabels=1, yticklabels=1, cmap='Blues')
        plt.ylabel('Organism isolates')
        plt.xlabel(x_axis_name)
        ax.figure.tight_layout()
        save_file = os.path.join(output_directory, image_name)
        plt.savefig(save_file, dpi=1200)
        plt.close()


def generate_bar_graph(image_name, x_axis_name, y_axis_name, default_dict_data_item, output_directory):
    """

    Generate a bar graph of the data. The default_dict_data_item is a dictionary with keys representing the
    resistance genes or plasmids. The lists of values are the isolates that contain the gene or plasmid.

    :param image_name: output file name
    :param default_dict_data_item: default dictionary with keys representing the resistance genes or plasmids
    :param output_directory: output directory
    """
    genes = []
    isolates = []
    combined_tuples = []
    for gene, isolate_list in default_dict_data_item.items():
        combined_tuples.append((gene, len(isolate_list)))
    sorted_tuples = sorted(combined_tuples, key=lambda x: x[1], reverse=True)
    for gene, isolate_count in sorted_tuples:
        genes.append(gene)
        isolates.append(isolate_count)
    sns.set(font_scale=0.8)
    sns.set_style("ticks")

    # Set the figure size
    plt.figure(figsize=(len(genes)*0.3, 10))
    sns.barplot(x=genes, y=isolates, palette="crest")
    # plt.bar(genes, isolates, width=0.3)
    plt.ylabel(y_axis_name)
    plt.xlabel(x_axis_name)
    plt.xticks(rotation='vertical')
    plt.tight_layout()
    save_file = os.path.join(output_directory, image_name)
    plt.savefig(save_file, dpi=600)
    plt.close()


def generate_pie_chart(image_name, default_dict, output_directory):
    """

    Generate a pie chart of the data.

    :param image_name: output file name
    :param default_dict: default dictionary with keys representing the type
    :param output_directory: output directory
    """
    type_list = []
    type_number = []
    for mlst_type, isolates in default_dict.items():
        type_list.append(mlst_type)
        type_number.append(len(isolates))

    sns.set(font_scale=0.5)
    plt.pie(type_number, labels=type_list, colors=sns.color_palette("Spectral", len(type_list)))
    plt.tight_layout()
    save_file = os.path.join(output_directory, image_name)
    plt.savefig(save_file, dpi=600)
    plt.close()

# def generate_polar_plot(image_name, default_dict, output_directory):
#     """
#
#     Generate a polar plot of the data.
#
#     :param image_name: output file name
#     :param default_dict: default dictionary with keys representing the type
#     :param output_directory: output directory
#     """
#     type_list = []
#     type_number = []
#     for mlst_type, isolates in default_dict.items():
#         type_list.append(mlst_type)
#         type_number.append(len(isolates))
#
#     N = len(type_list)
#     theta = np.linspace(0.0, 2 * np.pi, N, endpoint=False)
#     width = (2 * np.pi) / N
#
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='polar')
#
#     bars = ax.bar(theta, type_number, width=width, bottom=0.0)
#     ax.set_theta_zero_location('N')
#     ax.set_yticklabels([])
#
#     for r, bar in enumerate(bars):
#         angle = theta[r] + width/2
#         xpos = 0.5*bar.get_x()+0.5*bar.get_width()
#         ypos = bar.get_y()+bar.get_height()
#         ax.text(xpos, ypos, type_list[r], ha='center', va='bottom', rotation=angle)
#
#     plt.tight_layout()
#     save_file = os.path.join(output_directory, image_name)
#     plt.savefig(save_file, dpi=600)
#     plt.close()

# def generate_network_graph(image_name, default_dict_data_item, output_directory):
#     plasmid_list = []
#     genes_set = set()
#     labels = {}
#     for plasmid, genes in default_dict_data_item.items():
#         plasmid_list.append(plasmid)
#         labels[plasmid] = plasmid
#         for gene in genes:
#             labels[gene] = gene
#             genes_set.add(gene)
#
#     G = nx.Graph()
#     for plasmid in plasmid_list:
#         G.add_node(plasmid)
#     for gene in genes_set:
#         G.add_node(gene)
#
#     for plasmid, genes in default_dict_data_item.items():
#         for gene in genes:
#             G.add_edge(plasmid, gene)
#
#     pos = nx.circular_layout(G)
#     degrees = dict(G.degree())
#     node_sizes = [500*degrees[node] for node in G.nodes()]
#     nx.draw_networkx_nodes(G, pos, node_size=node_sizes)
#     nx.draw_networkx_edges(G, pos)
#     nx.draw_networkx_labels(G, pos, labels)
#     plt.axis('off')
#     save_file = os.path.join(output_directory, image_name)
#     plt.savefig(save_file)
#     plt.close()


def process_all_data(resfinder_results_directory,
                     plasmidfinder_result_directory,
                     mlst_result_directory,
                     quast_results_directory,
                     virulence_results_directory,
                     kaptive_results_directory):
    """
    This function populates data containers for secondary analysis of the results of the pipelines
    which output their data in the standard file format used by this pipeline.

    :param resfinder_results_directory: The main resfinder results directory
    :param plasmidfinder_result_directory: The main plasmidfinder results directory
    :param mlst_result_directory: The main mlst results directory
    :param quast_results_directory: The main quast results directory
    :return: Secondary analysis data containers
    """
    # Data dict with key the resistance gene and data: list of isolates with that gene
    processed_data_dict_isolate_res_genes = defaultdict(list)
    processed_data_dict_isolate_res_genes_org_separated = {}

    # Data dict with key the plasmid and data: list of resistance genes located in the same contig as that plasmid
    processed_data_dict_plasmid_resistance_genes = defaultdict(list)
    processed_data_dict_plasmid_resistance_genes_org_separated = {}

    # Data dict with key the plasmid and data: list of isolates with that plasmid
    processed_data_dict_isolate_plasmids = defaultdict(list)
    processed_data_dict_isolate_plasmids_org_separated = {}

    # Data dict with key the organism and the inner data a
    # defaultdict with key the type and data: list of isolates
    processed_data_dict_organism_capsule_type = {}

    # Data dict with key the isolate and data: list of virulence genes
    processed_data_dict_isolate_virulence = defaultdict(list)
    processed_data_dict_isolate_virulence_org_separated = {}
    processed_data_dict_isolate_virulence_simplified = defaultdict(list)
    processed_data_dict_isolate_virulence_org_separated_simplified = {}
    processed_data_dict_isolate_virulence_category = defaultdict(list)
    processed_data_dict_isolate_virulence_category_org_separated = {}

    # Data dict with key the isolate and data: list of contigs with a (1) resistance gene  (2) plasmid
    processed_data_dict_isolate_contig_resistance = {}
    processed_data_dict_isolate_contig_plasmid = {}

    general_data = defaultdict(list)
    mlst_data = defaultdict(list)
    mlst_all = {}
    isolate_genes = {}
    isolate_plasmids = {}
    isolate_quast = {}
    isolate_virulence_full = {}
    isolate_virulence_simplified = {}
    isolate_capsule = defaultdict(list)

    # Data dict with key the isolate and data: defaultdict with key the resistance gene and data:
    # list of contigs with that gene
    resistance_gene_contig_data = {}

    # Retrieve the datastructures from created files
    resfinder_data_structure = services.generate_result_info_dict(resfinder_results_directory)
    plasmidfinder_data_structure = services.generate_result_info_dict(plasmidfinder_result_directory)
    mlst_data_structure = services.generate_result_info_dict(mlst_result_directory)
    quast_data_structure = services.generate_result_info_dict(quast_results_directory)
    virulence_data_structure = services.generate_result_info_dict(virulence_results_directory)
    kaptive_data_structure = services.generate_result_info_dict(kaptive_results_directory)

    isolate_set = set()

    # Retrieve the mlst data
    for org, org_data_dicts in mlst_data_structure.items():
        org_data_dict = defaultdict(list)
        for org_data in org_data_dicts:
            isolate_id = org_data[1]
            isolate_set.add(isolate_id)
            file_path = os.path.join(org_data[0], 'data.json')
            processed_data = generate_mlst_datadicts(file_path)
            mlst_data[isolate_id].append(processed_data[0])
            mlst_data[isolate_id].append(processed_data[1])
            org_data_dict[processed_data[1]].append(isolate_id)
        mlst_all[org] = org_data_dict

    # Retrieve the resfinder data
    for org, org_data_dicts in resfinder_data_structure.items():
        org_data_dict = defaultdict(list)
        for org_data in org_data_dicts:
            isolate_id = org_data[1]
            isolate_set.add(isolate_id)
            file_path = os.path.join(org_data[0], 'ResFinder_results_tab.txt')
            processed_data = generate_resistance_datadicts(file_path)
            gene_contig_dict = defaultdict(list)
            all_genes = []
            for gene, isolate_data in processed_data[0].items():
                all_genes.append(gene)
                processed_data_dict_isolate_res_genes[gene].append(isolate_id)
                org_data_dict[gene].append(isolate_id)
                for isolation_data in isolate_data:
                    gene_contig_dict[isolation_data['Contig']].append(gene)

            isolate_genes[isolate_id] = all_genes
            resistance_gene_contig_data[isolate_id] = gene_contig_dict

        processed_data_dict_isolate_res_genes_org_separated[org] = org_data_dict

    # Retrieve the plasmidfinder data
    for org, org_data_dicts in plasmidfinder_data_structure.items():
        org_data_dict = defaultdict(list)
        org_data_dict_plasmid_isolates = defaultdict(list)
        for org_data in org_data_dicts:
            isolate_id = org_data[1]
            isolate_set.add(isolate_id)

            # Retrieve the plasmidfinder data from the file structure
            file_path = os.path.join(org_data[0], 'data.json')
            processed_data = generate_plasmid_datadicts(file_path)

            # Add the data to relevant data objects
            all_plasmids = []
            for plasmid, contig_data in processed_data[0].items():
                all_plasmids.append(plasmid)

                # Retrieve isolate data. A defaultdict with the key the contig
                # and data: list of resistance genes
                isolate_data = resistance_gene_contig_data[isolate_id]

                # Add the isolate to the defaultdict with key the plasmid and data: list of isolates
                # A total dict:
                processed_data_dict_isolate_plasmids[plasmid].append(isolate_id)
                # and a organism specific dict is populated:
                org_data_dict_plasmid_isolates[plasmid].append(isolate_id)

                # This iterates through the contig data and if the contig was the one used to
                # identify a resistance gene for this isolate it is added to the organism specific
                # org_data_dict and the total processed_data_dict_plasmid_resistance_genes
                for contig in contig_data:
                    if contig in isolate_data:
                        org_data_dict[plasmid].extend(isolate_data[contig])
                        processed_data_dict_plasmid_resistance_genes[plasmid].extend(isolate_data[contig])

            # The total data dict with key the isolate and data: list of plasmids
            isolate_plasmids[isolate_id] = all_plasmids

        # The organism separated data dicts are populated
        processed_data_dict_plasmid_resistance_genes_org_separated[org] = org_data_dict
        processed_data_dict_isolate_plasmids_org_separated[org] = org_data_dict_plasmid_isolates

    # Retrieve kaptive data
    for org, org_data_dicts in kaptive_data_structure.items():
        capsule_type_dict = {}
        for db_data in std_paths.kaptive_db_path[org]:
            capsule_defaultdict = defaultdict(list)
            for org_data in org_data_dicts:
                isolate_id = org_data[1]
                isolate_set.add(isolate_id)
                file_path = os.path.join(org_data[0], f'{db_data[0]}.json')
                processed_data = generate_kaptive_datadicts(file_path)
                if processed_data[0] != '':
                    capsule_defaultdict[processed_data[0]].append(isolate_id)
                    isolate_capsule[isolate_id].append(f'{processed_data[0]}: {processed_data[1]}')
            capsule_type_dict[db_data[0]] = capsule_defaultdict
        processed_data_dict_organism_capsule_type[org] = capsule_type_dict


    # Retrieve the quast data
    for org, org_data_dicts in quast_data_structure.items():
        for org_data in org_data_dicts:
            isolate_id = org_data[1]
            isolate_set.add(isolate_id)
            file_path = os.path.join(org_data[0], 'transposed_report.tsv')
            processed_data = generate_quast_datadicts(file_path)
            org_data_dict = {}
            for qc_key, qc_value in processed_data.items():
                org_data_dict[qc_key] = qc_value
            isolate_quast[isolate_id] = org_data_dict

    qc_values_report = services.quast_keys

    # Retrieve virulence data
    for org, org_data_dicts in virulence_data_structure.items():
        if org in services.vfdb_available:
            # Data dicts for organism separated data
            org_data_dict = defaultdict(list)
            org_data_dict_simple = defaultdict(list)
            org_data_dict_category = defaultdict(list)

            for org_data in org_data_dicts:
                isolate_id = org_data[1]
                isolate_set.add(isolate_id)

                # Retrieve virulence data from the file structure
                file_path = os.path.join(org_data[0], 'data.json')
                processed_data = virulence_data_handling.generate_virulence_datadicts(file_path, std_paths.virulence_info)

                # Add all virulence gene data to the full gene name data containers
                virulence_genes_all = []
                for vir_key, vir_value in processed_data[0].items():
                    virulence_genes_all.append(vir_key)
                    org_data_dict[vir_key].append(isolate_id)
                    processed_data_dict_isolate_virulence[vir_key].append(isolate_id)
                processed_data_dict_isolate_virulence_org_separated[org] = org_data_dict
                isolate_virulence_full[isolate_id] = virulence_genes_all

                # Add the virulence genes to the simplified name data containers
                virulence_genes_all_simple = []
                for virulence_simple in processed_data[2]:
                    virulence_genes_all_simple.append(virulence_simple)
                    org_data_dict_simple[virulence_simple].append(isolate_id)
                    processed_data_dict_isolate_virulence_simplified[virulence_simple].append(isolate_id)
                processed_data_dict_isolate_virulence_org_separated_simplified[org] = org_data_dict_simple
                isolate_virulence_simplified[isolate_id] = virulence_genes_all_simple

                # Add the virulence genes to categories
                for vir_cat_key, vir_cat_values in processed_data[3].items():
                    processed_data_dict_isolate_virulence_category[vir_cat_key].extend(vir_cat_values)
                    org_data_dict_category[vir_cat_key].extend(vir_cat_values)
                processed_data_dict_isolate_virulence_category_org_separated[org] = org_data_dict_category

    # Populate the general_data defaultdict which will be used to output the general Excel sheet
    for isolate_final in isolate_set:
        if isolate_final in mlst_data:
            general_data[isolate_final].append(mlst_data[isolate_final][0])
            general_data[isolate_final].append(mlst_data[isolate_final][1])
        else:
            general_data[isolate_final].append(services.org_name_key[isolate_final[:4]])
            general_data[isolate_final].append('No schema available')

        if isolate_final in isolate_capsule:
            all_capsules = isolate_capsule[isolate_final]
            if len(all_capsules) > 1:
                capsule_string = ', '.join(all_capsules)
                general_data[isolate_final].append(capsule_string)
            elif len(isolate_capsule[isolate_final]) == 1:
                general_data[isolate_final].append(str(all_capsules[0]))
            else:
                general_data[isolate_final].append('None')

        if isolate_final in isolate_genes:
            all_genes = isolate_genes[isolate_final]
            if len(all_genes) > 1:
                gene_string = ', '.join(all_genes)
                general_data[isolate_final].append(gene_string)
            elif len(isolate_genes[isolate_final]) == 1:
                general_data[isolate_final].append(str(all_genes[0]))
            else:
                general_data[isolate_final].append('None')
        else:
            general_data[isolate_final].append('None')

        if isolate_final in isolate_plasmids:
            all_plasmids = isolate_plasmids[isolate_final]
            if len(all_plasmids) > 1:
                plasmids_string = ', '.join(all_plasmids)
                general_data[isolate_final].append(plasmids_string)
            elif len(all_plasmids) == 1:
                general_data[isolate_final].append(str(all_plasmids[0]))
            else:
                general_data[isolate_final].append('None')
        else:
            general_data[isolate_final].append('None')

        if isolate_final in isolate_virulence_full:
            all_virulence = isolate_virulence_full[isolate_final]
            if len(all_virulence) > 1:
                virulence_string = ', '.join(all_virulence)
                general_data[isolate_final].append(virulence_string)
            elif len(all_virulence) == 1:
                general_data[isolate_final].append(str(all_virulence[0]))
            else:
                general_data[isolate_final].append('None')
        else:
            general_data[isolate_final].append('None')

        if isolate_final in isolate_quast:
            for qc_val in qc_values_report:
                if qc_val in isolate_quast[isolate_final]:
                    general_data[isolate_final].append(isolate_quast[isolate_final][qc_val])
                else:
                    general_data[isolate_final].append('Not evaluated')
        else:
            for _ in qc_values_report:
                general_data[isolate_final].append('Not evaluated')

    return (processed_data_dict_isolate_res_genes,
            processed_data_dict_isolate_res_genes_org_separated,
            processed_data_dict_plasmid_resistance_genes,
            processed_data_dict_plasmid_resistance_genes_org_separated,
            general_data,
            isolate_genes,
            isolate_plasmids,
            qc_values_report,
            mlst_all,
            processed_data_dict_isolate_plasmids,
            processed_data_dict_isolate_plasmids_org_separated,
            processed_data_dict_isolate_virulence,
            processed_data_dict_isolate_virulence_org_separated,
            isolate_virulence_full,
            isolate_virulence_simplified,
            processed_data_dict_isolate_virulence_category,
            processed_data_dict_isolate_virulence_category_org_separated,
            processed_data_dict_isolate_virulence_simplified,
            processed_data_dict_isolate_virulence_org_separated_simplified,
            processed_data_dict_organism_capsule_type)


def generate_figures_data(output_directory, resfinder_results_directory, plasmidfinder_result_directory,
                          mlst_result_directory, quast_result_directory, virulence_result_directory,
                          kaptive_result_directory):

    # Process the data for display
    processed_data = process_all_data(resfinder_results_directory,
                                      plasmidfinder_result_directory,
                                      mlst_result_directory,
                                      quast_result_directory,
                                      virulence_result_directory,
                                      kaptive_result_directory)

    # Generate the general and organism specific resistance gene heatmaps
    generate_heat_map('All_organisms_heatmap.png', 'Resistance genes',
                      processed_data[0], output_directory)
    for org, data in processed_data[1].items():
        generate_heat_map(f'{org}_heatmap.png', 'Resistance genes',
                          data, output_directory)

    # Generate the general and organism specific plasmid resistance gene heatmaps
    generate_heat_map('All_organisms_plasmid_heatmap.png', 'Plasmids',
                      processed_data[9], output_directory)
    for org, data in processed_data[10].items():
        generate_heat_map(f'{org}_plasmid_heatmap.png', 'Plasmids',
                          data, output_directory)

    # Generate bar graph for resistance genes
    generate_bar_graph('All_organisms_resistance_gene_bar_graph.png', 'Resistance genes', 'Number of isolates',
                       processed_data[0], output_directory)
    for org, data in processed_data[1].items():
        generate_bar_graph(f'{org}_resistance_gene_bar_graph.png', 'Resistance genes', 'Number of isolates',
                           data, output_directory)

    # Generate bar graph for plasmids
    generate_bar_graph('All_organisms_plasmid_bar_graph.png', 'Plasmids', 'Number of isolates',
                       processed_data[9], output_directory)
    for org, data in processed_data[10].items():
        generate_bar_graph(f'{org}_plasmid_bar_graph.png', 'Plasmids', 'Number of isolates',
                           data, output_directory)

    # Generate bar graph for virulence genes
    generate_bar_graph('All_organisms_virulence_gene_bar_graph.png', 'Virulence genes', 'Number of genes detected',
                       processed_data[17], output_directory)
    for org, data in processed_data[18].items():
        generate_bar_graph(f'{org}_virulence_gene_bar_graph.png', 'Virulence genes', 'Number of genes detected',
                           data, output_directory)

    # Generate pie chart for virulence genes
    pie_chart_dict_all = {}
    for category, virulence_genes in processed_data[15].items():
        virulence_set = set()
        for gene in virulence_genes:
            virulence_set.add(gene)
        pie_chart_dict_all[category] = virulence_set
    generate_pie_chart('All_organisms_virulence_gene_category_pie_chart.png',
                       pie_chart_dict_all, output_directory)
    for org, data in processed_data[16].items():
        pie_chart_dict = {}
        for category, virulence_genes in data.items():
            virulence_set = set()
            for gene in virulence_genes:
                virulence_set.add(gene)
            pie_chart_dict[category] = virulence_set
        generate_pie_chart(f'{org}_virulence_gene_category_pie_chart.png',
                           pie_chart_dict, output_directory)

    # Generate pie chart for MLST for each organism
    for org, data in processed_data[8].items():
        generate_pie_chart(f'{org}_mlst_pie_chart.png', data, output_directory)

    # Generate pie chart for capsule type for each organism
    for org, data in processed_data[19].items():
        for capsule_type, isolates_def_dict in data.items():
            if len(isolates_def_dict.items()) > 1:
                generate_pie_chart(f'{org}_{capsule_type}_capsule_type_pie_chart.png', isolates_def_dict,
                                   output_directory)

    # Generate excel file with all resistance genes
    resistance_gene_data = processed_data[1]
    for org, org_data in resistance_gene_data.items():
        output_dict = defaultdict(list)
        resistance_genes = []
        isolate_set = set()
        for resistance_gene, isolate_data in org_data.items():
            resistance_genes.append(resistance_gene)
            for isolate in isolate_data:
                isolate_set.add(isolate)
        resistance_genes_sorted = sorted(resistance_genes)
        for isolate in isolate_set:
            for resistance_gene in resistance_genes_sorted:
                if isolate in org_data[resistance_gene]:
                    output_dict[isolate].append('1')
                else:
                    output_dict[isolate].append('0')
        df = pd.DataFrame(output_dict, index=resistance_genes_sorted).transpose()
        writer = pd.ExcelWriter(os.path.join(output_directory, f'{org}: Resistance_raw_data.xlsx'))
        df.to_excel(writer)
        writer.close()

    # Generate excel file with all plasmid data
    plasmid_data = processed_data[10]
    for org, org_data in plasmid_data.items():
        output_dict = defaultdict(list)
        plasmids = []
        isolate_set = set()
        for plasmid, isolate_data in org_data.items():
            plasmids.append(plasmid)
            for isolate in isolate_data:
                isolate_set.add(isolate)
        plasmids_sorted = sorted(plasmids)
        for isolate in isolate_set:
            for plasmid_detect in plasmids_sorted:
                if isolate in org_data[plasmid_detect]:
                    output_dict[isolate].append('1')
                else:
                    output_dict[isolate].append('0')
        df = pd.DataFrame(output_dict, index=plasmids_sorted).transpose()
        writer = pd.ExcelWriter(os.path.join(output_directory, f'{org}: Plasmids_raw_data.xlsx'))
        df.to_excel(writer)
        writer.close()

    # Generate excel file with all virulence data
    virulence_data = processed_data[12]
    for org, org_data in virulence_data.items():
        output_dict = defaultdict(list)
        isolates = []
        virulence_set = set()
        for isolate, virulence_data in org_data.items():
            isolates.append(isolate)
            for virulence_gene in virulence_data:
                virulence_set.add(virulence_gene)
        for gene in virulence_set:
            for isolate in isolates:
                if gene in org_data[isolate]:
                    output_dict[gene].append('1')
                else:
                    output_dict[gene].append('0')
        df = pd.DataFrame(output_dict, index=isolates).transpose()
        writer = pd.ExcelWriter(os.path.join(output_directory, f'{org}: Virulence_raw_data.xlsx'))
        df.to_excel(writer)
        writer.close()


    # Output results to an Excel file
    output_table_headings = ['Organism', 'MLST Sequence Type', 'Capsule Type', 'Resistance genes', 'Plasmids',
                             'Virulence genes']
    output_table_headings.extend(processed_data[7])
    df = pd.DataFrame(processed_data[4], index=output_table_headings).transpose()
    writer = pd.ExcelWriter(os.path.join(output_directory, 'General_data_summary.xlsx'))
    df.to_excel(writer)
    writer.close()

    # generate
    # generate_network_graph('All_organisms_network.png', processed_data[2], output_directory)
    # for org, data in processed_data[3].items():
        # generate_network_graph(f'{org}_network.png', data, output_directory)


# print(generate_mlst_datadicts('test/mlst_results/KLEB/KLEB-CRE-GSH-0001/data.json'))



