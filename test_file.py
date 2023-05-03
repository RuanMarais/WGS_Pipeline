import json
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

print(generate_kaptive_datadicts('test/kaptive_results/KLEB/KLEB-CRE-GSH-0001/K_primary.json'))