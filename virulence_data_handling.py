from collections import defaultdict
import json
import pandas as pd


def generate_virulence_info(data_file):
    # read the Excel file into a pandas DataFrame, skipping the first row
    df = pd.read_excel(data_file, skiprows=1)
    df.set_index(df.columns[0], inplace=True)
    return df


def virulence_name_and_id(fasta_header):
    bracket_index = fasta_header.find('[')
    full_details = fasta_header[:bracket_index-1]
    post_bracket = fasta_header[bracket_index+1:]
    vf_start = post_bracket.find('VF')
    VF_id = post_bracket[vf_start:vf_start+6]
    vf_short = post_bracket[:vf_start-2]
    vf_fx_start = post_bracket.find('_-_')
    vf_fx_end = post_bracket.find('_(VFC')
    vf_fx = post_bracket[vf_fx_start+3:vf_fx_end]
    return full_details, VF_id, vf_short, vf_fx


def generate_virulence_datadicts(data_file, virulence_info):
    output_dict_virulence_sorted = defaultdict(list)
    output_dict_contig_sorted = defaultdict(list)
    simplified_virulence = []
    output_category_dict = defaultdict(list)

    try:
        with open(data_file, 'r') as f:
            json_string = f.read()
            data = json.loads(json_string)
            data_items = data['mydbfinder']['results']['Vfdb_all']['VFDB_all']
            for key, val in data_items.items():
                virulence_data = virulence_name_and_id(val['fasta_header'])
                vf_category = virulence_data[3]
                output_dict_virulence_sorted[virulence_data[0]].append((virulence_data, vf_category, val['contig_name']))
                simplified_virulence.append(f'{virulence_data[2]}: {vf_category}')
                output_dict_contig_sorted[val['contig_name']].append(val['fasta_header'])
                output_category_dict[vf_category].append(f'{virulence_data[2]}: {vf_category}')
    except:
        with open('test/error_log', 'a+') as f:
            f.seek(0, 2)
            f.write(f'Generate virulence datadicts failed for: {data_file}\n')
    return output_dict_virulence_sorted, output_dict_contig_sorted, simplified_virulence, output_category_dict
