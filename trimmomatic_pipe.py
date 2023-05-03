import os


def generate_trimmomatic_commands(paired_file_list, id_length, trim_results_directory, trimmomatic_path, threads,
                                  source_folder_path, adaptors_path):
    output = []
    for paired_files in paired_file_list:
        id_file = paired_files[0][:id_length + 1]
        result_paired_1 = os.path.join(trim_results_directory, id_file + '_paired_1.fastq')
        result_paired_2 = os.path.join(trim_results_directory, id_file + '_paired_2.fastq')
        result_unpaired_1 = os.path.join(trim_results_directory, id_file + '_unpaired_1.fastq')
        result_unpaired_2 = os.path.join(trim_results_directory, id_file + '_unpaired_2.fastq')
        if trimmomatic_path is None:
            trim_command = ['trimmomatic', 'PE', '-threads', str(threads), source_folder_path + '/' + paired_files[0],
                             source_folder_path + '/' + paired_files[1], result_paired_1, result_unpaired_1,
                             result_paired_2, result_unpaired_2, f'ILLUMINACLIP:{adaptors_path}:2:30:10',
                             'LEADING:20', 'TRAILING:20', 'SLIDINGWINDOW:4:20', 'MINLEN:50']
        else:
            trim_command = ['java', '-jar', trimmomatic_path, 'PE', '--threads', str(threads),
                             source_folder_path + '/' + paired_files[0], source_folder_path + '/' + paired_files[1],
                             result_paired_1, result_unpaired_1, result_paired_2, result_unpaired_2,
                             f'ILLUMINACLIP:{adaptors_path}:2:30:10', 'LEADING:20', 'TRAILING:20',
                             'SLIDINGWINDOW:4:20', 'MINLEN:50']
        output.append((id_file, trim_command))
    return output
