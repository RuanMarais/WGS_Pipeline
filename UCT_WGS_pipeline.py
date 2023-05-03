import subprocess
import BRIG_file_creation
import fasttree_pipe
import gen_find
import pilon_pipe
import reference_selection
import services
import argparse
import os
from collections import defaultdict
import spades_pipe
import trimmomatic_pipe
import resfinder_pipe
import plasmid_finder_pipe
import data_handling
import mlst_pipe
import quast_pipe
import prokka_pipe
import std_paths
import parsnp_pipe
import kmerfinder_pipe
import fastqc_pipe
import ragtag_pipe
import kaptive_pipe
import hts_functions
import mgefinder_pipe
import shutil
import resfinder_mgeresults
import isescan_mgeresults
import roary_pipe
import roary_plot_pipe
import ncbi_refseq_download_pipe
import poppunk_pipe
import kmeans_clustering
import statistics
import reference_evaluation


# The input file for the pipeline
parser = argparse.ArgumentParser(description='This pipeline analyses WGS bacterial data. Created by GJK Marais.')
parser.add_argument('--input_folder', required=True, help='The file path that contains Illumina paired-end reads')
parser.add_argument('--output_folder', required=True, help='Result folders will be output to this folder')
parser.add_argument('--threads', default=1, help='The number of threads to use, default is 1')
parser.add_argument('--trimmomatic', default=None, type=int, help='Path to trimmomatic JAR file. Default is bin')
parser.add_argument('--trim_adaptors_path', default='TruSeq3-PE.fa', help='Path to the adaptors to remove')
parser.add_argument('--identifier_length', type=int, required=True,
                    help='The string length that identifies the sequencing read')
parser.add_argument('--organism', default=None, help='The organism >Genus species< of provided raw data. '
                                                     'If none provided, it will be detected from the data')
parser.add_argument('--file_size_deviation', default=0.2, type=float,
                    help='The allowable file_size deviation for refseq downloaded sequences')
parser.add_argument('--download', action='store_true',
                    help='If selected, the pipeline will download the reference genomes from Refseq')
parser.add_argument('--path_to_refseq_genomes', default=None, help='Path to Refseq genomes to use select a reference')
args = parser.parse_args()

# Set variables for working directory, id_length and threads
input_folder = args.input_folder
folder_path = args.output_folder
id_length = args.identifier_length
threads = args.threads
organism_name = args.organism
file_deviation = args.file_size_deviation
download = args.download
refseq_genomes = args.path_to_refseq_genomes

# Initialise the error file data item
error_log = []
error_file = 'error_log'
if os.path.isfile(error_file):
    os.remove(os.path.join(folder_path, error_file))

# Initialise the log file data item
log_file = 'log_file'
if os.path.isfile(log_file):
    os.remove(os.path.join(folder_path, log_file))

# Generate a list with all valid filenames
file_names_all = services.filename_list_generate('001.fastq.gz', input_folder)

# Separate filenames by id and into a list of pair tuples
paired_sorted_list = services.paired_read_list_generate(id_length, 'R1_001.fastq.gz',
                                                        'R2_001.fastq.gz', file_names_all)

# Variable for the list of trimming commands
trim_commands_list = []

# Create the results directories
trim_results_directory = os.path.join(folder_path, 'trim_results')
fastqc_results_directory = os.path.join(folder_path, 'fastqc_results')
spades_results_directory = os.path.join(folder_path, 'spades_results')
pilon_results_directory = os.path.join(folder_path, 'pilon_results')
quast_results_directory = os.path.join(folder_path, 'quast_results')
resfinder_results_directory = os.path.join(folder_path, 'resfinder_results')
plasmidfinder_results_directory = os.path.join(folder_path, 'plasmidfinder_results')
figure_results_directory = os.path.join(folder_path, 'figures')
mlst_results_directory = os.path.join(folder_path, 'mlst_results')
prokka_results_directory = os.path.join(folder_path, 'prokka_results')
virulence_results_directory = os.path.join(folder_path, 'virulence_results')
parsnp_evaluation_directory = os.path.join(folder_path, 'parsnp_evaluation')
kmerfinder_results_directory = os.path.join(folder_path, 'kmerfinder_results')
ragtag_results_directory = os.path.join(folder_path, 'ragtag_results')
kaptive_results_directory = os.path.join(folder_path, 'kaptive_results')
deduplicated_results_directory = os.path.join(folder_path, 'deduplicated_results')
mge_directory = os.path.join(folder_path, 'mge_directory')
isescan_results_directory = os.path.join(folder_path, 'isescan_results')
mge_resfinder_results_directory = os.path.join(folder_path, 'mge_resfinder_results')
roary_results_directory = os.path.join(folder_path, 'roary_results')
fasttree_results_directory = os.path.join(folder_path, 'fasttree_results')
if refseq_genomes is None:
    refseq_genomes = os.path.join(folder_path, 'refseq_genomes')
    if not os.path.exists(refseq_genomes):
        os.makedirs(refseq_genomes)
candidate_reference_genomes = os.path.join(folder_path, 'candidate_reference_genomes')
candidate_reference_genome_alignment = os.path.join(folder_path, 'candidate_reference_genome_alignment')
candidate_reference_genome_alignment_output = os.path.join(folder_path, 'candidate_reference_genome_alignment_output')
poppunk_results_reference = os.path.join(folder_path, 'poppunk_results_reference')
selected_reference_genome = os.path.join(folder_path, 'selected_reference_genome')
brig_evaluation_directory = os.path.join(folder_path, 'brig_evaluation')

directories_to_create = [trim_results_directory, fastqc_results_directory, spades_results_directory,
                         pilon_results_directory, quast_results_directory, resfinder_results_directory,
                         plasmidfinder_results_directory, figure_results_directory, mlst_results_directory,
                         prokka_results_directory, virulence_results_directory, parsnp_evaluation_directory,
                         kmerfinder_results_directory, ragtag_results_directory, kaptive_results_directory,
                         deduplicated_results_directory, mge_directory, isescan_results_directory,
                         mge_resfinder_results_directory, roary_results_directory, fasttree_results_directory,
                         candidate_reference_genomes, candidate_reference_genome_alignment,
                         candidate_reference_genome_alignment_output, poppunk_results_reference,
                         selected_reference_genome, brig_evaluation_directory]

for directory in directories_to_create:
    if not os.path.exists(directory):
        os.makedirs(directory)

# Generate deduplicated reads
dedup_commands = hts_functions.generate_deduplicator_commands(paired_sorted_list, id_length, input_folder,
                                                              deduplicated_results_directory)

# Run the deduplicator commands with error handling
for command_item in dedup_commands:
    try:
        subprocess.run(command_item, check=True)
    except subprocess.CalledProcessError as e:
        error_log.append(f'hts_SuperDeduper error for command: {command_item} with error code: {e.returncode}')

# Get all files in deduplicated folder
file_names_to_trim = services.filename_list_generate('fastq.gz', deduplicated_results_directory)
paired_list_for_trim = services.paired_read_list_generate(id_length, '_R1.fastq.gz', '_R2.fastq.gz', file_names_to_trim)

# Create the trimmomatic commands: Read QC
trim_commands_list = trimmomatic_pipe.generate_trimmomatic_commands(paired_list_for_trim,
                                                                    id_length,
                                                                    trim_results_directory,
                                                                    args.trimmomatic,
                                                                    threads,
                                                                    deduplicated_results_directory,
                                                                    args.trim_adaptors_path)

# Run the trimmomatic commands with error handling
for trim_command_item in trim_commands_list:
    try:
        subprocess.run(trim_command_item[1], check=True)
    except subprocess.CalledProcessError as e:
        error_log.append(f'trimmomatic error for specimen: {trim_command_item[0]} with error code: {e.returncode}')

# Trimmed file list
file_names_trimmed = services.filename_list_generate('.fastq', trim_results_directory)

# Separate filenames by id and into a list of pair tuples
paired_trimmed_list = services.paired_read_list_generate(id_length, '_paired_1.fastq',
                                                         '_paired_2.fastq', file_names_trimmed)

# Generate FastQC reports
fastqc_commands = fastqc_pipe.generate_fastqc_commands(paired_trimmed_list, id_length, fastqc_results_directory,
                                                       threads, trim_results_directory)
# Run the FastQC commands with error handling
for fastqc_command_item in fastqc_commands:
    try:
        subprocess.run(fastqc_command_item[1], check=True)
    except subprocess.CalledProcessError as e:
        error_log.append(f'SPAdes error for specimen: {fastqc_command_item[0]} with error code: {e.returncode}')

# SPAdes commands: Assembly
spades_commands = spades_pipe.generate_spades_commands(paired_trimmed_list, id_length, threads,
                                                       spades_results_directory,
                                                       trim_results_directory)

# Run the SPAdes commands with error handling
for spades_command_item in spades_commands:
    try:
        subprocess.run(spades_command_item[1], check=True)
    except subprocess.CalledProcessError as e:
        error_log.append(f'SPAdes error for specimen: {spades_command_item[0]} with error code: {e.returncode}')

# Collect SPAdes generated contigs for further analysis
assemblies_to_evaluate = []
folders_spades = os.listdir(spades_results_directory)
directory_paths_spades = [(os.path.join(spades_results_directory, folder), folder)
                          for folder in folders_spades if os.path.isdir(os.path.join(spades_results_directory, folder))]
for directory in directory_paths_spades:
    contigs_detected = os.path.isfile(os.path.join(directory[0], 'contigs.fasta'))
    scaffold_detected = os.path.isfile(os.path.join(directory[0], 'scaffolds.fasta'))
    if contigs_detected and scaffold_detected:
        id_sample = directory[1]
        contigs_file = 'contigs.fasta'
        scaffold_file = 'scaffolds.fasta'
        assemblies_to_evaluate.append((id_sample, os.path.join(directory[0], scaffold_file),
                                       os.path.join(directory[0], contigs_file)))

# Pilon commands: Polishing
pilon_commands = pilon_pipe.generate_pilon_commands(assemblies_to_evaluate,
                                                    pilon_results_directory,
                                                    threads,
                                                    folder_path)

# Run the Pilon commands with error handling
for command_all in pilon_commands:
    for bowtie_commands in command_all[1][0]:
        try:
            subprocess.run(bowtie_commands, check=True)
        except subprocess.CalledProcessError as e:
            error_log.append(f'Bowtie2 error for specimen: {command_all[0]} with error code: {e.returncode}')
    try:
        subprocess.run(command_all[1][1], check=True)
    except subprocess.CalledProcessError as e:
        error_log.append(f'Pilon error for specimen: {command_all[0]} with error code: {e.returncode}')

assemblies_polished = []
folders_pilon = os.listdir(pilon_results_directory)
directory_paths_pilon = [(os.path.join(pilon_results_directory, folder), folder)
                         for folder in folders_pilon if os.path.isdir(os.path.join(pilon_results_directory, folder))]
for directory in directory_paths_pilon:
    file_detected = os.path.isfile(os.path.join(directory[0], directory[1] + '.fasta'))
    if file_detected:
        id_sample = directory[1]
        polished_file = id_sample + '.fasta'
        assemblies_polished.append((id_sample, os.path.join(directory[0], polished_file)))

if organism_name is None:
    # Generate kmerfinder commands
    kmerfinder_commands = kmerfinder_pipe.generate_kmerfinder_commands(assemblies_polished,
                                                                       kmerfinder_results_directory)

    # Run the kmerfinder commands with error handling to establish species
    for command in kmerfinder_commands:
        try:
            subprocess.run(command[1], check=True)
        except subprocess.CalledProcessError as e:
            error_log.append(f'{command[0]}: Error for command: {command[1]} with error code: {e.returncode}')

    # Generate organism dictionary
    organism_id_dictionary = {}
else:
    # Generate organism dictionary
    organism_id_dictionary = {}
    for assembly in assemblies_polished:
        organism_id_dictionary[assembly[0]] = organism_name

for assembly_org in organism_id_dictionary.items():
    # shutil copy the assembly to the candidate_reference_genome_alignment folder with subfolders for organisms
    for item in assemblies_polished:
        if assembly_org[0] == item[0]:
            organism_folder = os.path.join(candidate_reference_genome_alignment, assembly_org[1])
            if not os.path.exists(organism_folder):
                os.makedirs(organism_folder)
            shutil.copy(item[1], os.path.join(organism_folder, item[0] + '.fasta'))

if download:
    # Retrieve refseq sequences
    orgs_list_to_download_genomes = set()

    for key, org in organism_id_dictionary.items():
        orgs_list_to_download_genomes.add(org)

    refseq_commands = ncbi_refseq_download_pipe.download_genomes(refseq_genomes, orgs_list_to_download_genomes)

    # Run the refseq commands with error handling
    for command in refseq_commands:
        try:
            subprocess.run(command, check=True, shell=True)
        except subprocess.CalledProcessError as e:
            error_log.append(f'Error for command: {command} with error code: {e.returncode}')

# Retrieve reference genome folders
refseq_genome_folder_list = []

refseq_folder_names = os.listdir(refseq_genomes)
directory_paths_refseq = [(folder, os.path.join(refseq_genomes, folder))
                          for folder in refseq_folder_names if os.path.isdir(os.path.join(refseq_genomes, folder))]

# extract downloaded genomes
for directory in directory_paths_refseq:
    files = os.listdir(directory[1])
    for file in files:
        if file.endswith('.gz'):
            type_output_directory = os.path.join(poppunk_results_reference, directory[0])
            if not os.path.exists(type_output_directory):
                os.makedirs(type_output_directory)
            services.decompress_gzip_file(os.path.join(directory[1], file),
                                          os.path.join(type_output_directory, file[:-3]),
                                          error_log)

reference_genome_folder_list = []
reference_folder_names = os.listdir(poppunk_results_reference)
directory_paths_poppunk = [(folder, os.path.join(poppunk_results_reference, folder))
                           for folder in reference_folder_names
                           if os.path.isdir(os.path.join(poppunk_results_reference, folder))]

# Generate the filename text file for poppunk
# Files with irregular sizes are excluded
poppunk_filename_genomes = "filename_list.txt"
poppunk_genome_name_key = {}
for directory in directory_paths_poppunk:
    # get a list of all the filenames in the folder
    filenames = os.listdir(directory[1])
    filepaths = [os.path.join(directory[1], file)
                 for file in filenames if os.path.isfile(os.path.join(directory[1], file))]
    for file_to_format in filenames:
        poppunk_genome_name_key[services.format_filename_to_dashes(file_to_format)] = os.path.join(directory[1],
                                                                                                   file_to_format)
    file_sizes = [os.path.getsize(f) for f in filepaths]
    mean_file_size = statistics.mean(file_sizes)
    allowed_deviation = file_deviation * mean_file_size

    # create a new text file
    with open(os.path.join(directory[1], poppunk_filename_genomes), "w") as file:
        # loop through the filenames and write them to the text file
        for filename in filenames:
            filepath = os.path.join(directory[1], filename)
            if abs(os.path.getsize(filepath) - mean_file_size) < allowed_deviation:
                file.write(filename + "\t" + filepath + "\n")

# Generate popPUNK commands
poppunk_reference_eval_commands = poppunk_pipe.generate_poppunk_commands(directory_paths_poppunk,
                                                                         poppunk_filename_genomes,
                                                                         threads)

# Run the poppunk commands with error handling
for command in poppunk_reference_eval_commands:
    try:
        subprocess.run(command[0], check=True)
    except subprocess.CalledProcessError as e:
        error_log.append(f'Error for command: {command[0]} with error code: {e.returncode}')
    try:
        subprocess.run(command[1], check=True)
    except subprocess.CalledProcessError as e:
        error_log.append(f'Error for command: {command[1]} with error code: {e.returncode}')

# Run the kmeans clustering algorithm on all the genomes
reference_genome_defaultdict = defaultdict(list)
for directory_data in directory_paths_poppunk:
    distances = os.path.join(os.path.join(directory_data[1], 'poppunk_clusters'), 'poppunk_clusters.dists.out')
    genomes = kmeans_clustering.return_reference_genomes(distances, directory_data[1])
    org_output_directory = os.path.join(candidate_reference_genomes, directory_data[0])
    if not os.path.exists(org_output_directory):
        os.makedirs(org_output_directory)
    for genome in genomes:
        new_path_genome = os.path.join(org_output_directory, genome)
        shutil.copy(poppunk_genome_name_key[genome], new_path_genome)
        reference_genome_defaultdict[directory_data[0]].append((org_output_directory, new_path_genome, genome))


# Generate commands to index candidate references
bwa_index_candidates = reference_selection.generate_bwa_index(reference_genome_defaultdict)

# Run the bwa index commands with error handling
for command in bwa_index_candidates:
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        error_log.append(f'Error for command: {command} with error code: {e.returncode}')

# Generate file structure for alignment of assemblies to candidate references
assemblies_and_references = {}
folders_assemblies_to_align = os.listdir(candidate_reference_genome_alignment)
directory_paths_orgs = [(os.path.join(candidate_reference_genome_alignment, folder), folder)
                        for folder in folders_assemblies_to_align
                        if os.path.isdir(os.path.join(candidate_reference_genome_alignment, folder))]
for org_directory in directory_paths_orgs:
    path_to_genomes_to_align = os.path.join(candidate_reference_genome_alignment, org_directory[1])
    files_to_align = os.listdir(path_to_genomes_to_align)
    files_paths = [(file, os.path.join(path_to_genomes_to_align, file))
                   for file in files_to_align
                   if os.path.isfile(os.path.join(path_to_genomes_to_align, file))]
    references = reference_genome_defaultdict[org_directory[1]]
    assemblies_and_references[org_directory[1]] = [references, files_paths]

# Generate commands to align assemblies to candidate references
bwa_align_candidates = reference_selection.align_scaffolds_bwa(assemblies_and_references, threads,
                                                               candidate_reference_genome_alignment_output)
reference_selection.generate_consensus_sequences(bwa_align_candidates, error_log)

# Evaluate candidate references based on average identity
org_reference_folders = os.listdir(candidate_reference_genome_alignment_output)
directory_paths_orgs = [(os.path.join(candidate_reference_genome_alignment_output, folder), folder)
                        for folder in org_reference_folders
                        if os.path.isdir(os.path.join(candidate_reference_genome_alignment_output, folder))]

# This code evaluates the candidate reference genome average alignments and creates a dictionary for downstream use
selected_reference_path_dictionary = {}
for org_directory in directory_paths_orgs:
    output_reference_org_directory = os.path.join(selected_reference_genome, org_directory[1])
    if not os.path.exists(output_reference_org_directory):
        os.makedirs(output_reference_org_directory)
    candidate_references = os.listdir(org_directory[0])
    directories_to_evaluate = [(os.path.join(org_directory[0], folder), folder)
                               for folder in candidate_references
                               if os.path.isdir(os.path.join(org_directory[0], folder))]
    average_identity_count = 0
    selected_reference = None
    for genome_directory in directories_to_evaluate:
        average_identity = reference_evaluation.evaluate_average_identity(genome_directory[0])
        if average_identity > average_identity_count or selected_reference is None:
            average_identity_count = average_identity
            selected_reference = genome_directory[1]

    selected_reference_path = os.path.join(os.path.join(candidate_reference_genomes, org_directory[1]),
                                           selected_reference)
    reference_copy_path = os.path.join(output_reference_org_directory, selected_reference)
    shutil.copy(selected_reference_path, reference_copy_path)
    selected_reference_path_dictionary[org_directory[1]] = (reference_copy_path, selected_reference)


# Generate bwa alignments for mgefinder
mge_commands = mgefinder_pipe.generate_bwa_commands(paired_trimmed_list,
                                                    id_length,
                                                    trim_results_directory,
                                                    mge_directory,
                                                    threads,
                                                    organism_id_dictionary,
                                                    selected_reference_path_dictionary)

# Run mgefinder commands
for index_command in mge_commands[0]:
    try:
        subprocess.run(index_command, check=True)
    except subprocess.CalledProcessError as e:
        error_log.append(f'Database index for BWA: {index_command} with error code: {e.returncode}')

for align_command in mge_commands[1]:
    try:
        subprocess.run(align_command[0], check=True)
    except subprocess.CalledProcessError as e:
        error_log.append(f'Alignment for BWA: {align_command} with error code: {e.returncode}')
    try:
        subprocess.run(align_command[1], check=True)
    except subprocess.CalledProcessError as e:
        error_log.append(f'Sorting for BWA: {align_command} with error code: {e.returncode}')

# Sort contigs to source folder
for sample in assemblies_to_evaluate:
    org = organism_id_dictionary[sample[0]]
    shutil.copy(sample[2], os.path.join(os.path.join(os.path.join(mge_directory, org), '00.assembly'),
                                        f'{sample[0]}.fna'))

# Run mgefinder main commands
for mge_command in mge_commands[2]:
    try:
        subprocess.run(mge_command, check=True)
    except subprocess.CalledProcessError as e:
        error_log.append(f'Mgefinder error for specimen: {mge_command} with error code: {e.returncode}')

# Retrieve mgefinder results
mge_finder_results = {}
org_folders = os.listdir(mge_directory)
for org in org_folders:
    all_seq = os.path.join(os.path.join(os.path.join(os.path.join(mge_directory, org), '03.results'), org),
                           f'04.makefasta.{org}.all_seqs.fna')
    repr_seq = os.path.join(os.path.join(os.path.join(os.path.join(mge_directory, org), '03.results'), org),
                            f'04.makefasta.{org}.repr_seqs.fna')
    if os.path.isfile(all_seq) and os.path.isfile(repr_seq):
        mge_finder_results[org] = (all_seq, repr_seq)

# Create mgefinder results analysis commands: resfinder and isescan
resfinder_mge_commands = resfinder_mgeresults.generate_resfinder_commands(mge_finder_results,
                                                                          mge_resfinder_results_directory)
isescan_mge_commands = isescan_mgeresults.generate_isescan_commands(mge_finder_results,
                                                                    isescan_results_directory,
                                                                    threads)

# Run resfinder and isescan on mgefinder results
for command in resfinder_mge_commands:
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        error_log.append(f'Resfinder error for MGEfinder results: {command} with error code: {e.returncode}')

for command in isescan_mge_commands:
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        error_log.append(f'ISEScan error for MGEfinder results: {command} with error code: {e.returncode}')

# Generate org separated assemblies dictionary for ragtag
org_separated_polished_assemblies = defaultdict(list)
for assembly in assemblies_polished:
    org = organism_id_dictionary[assembly[0]]
    org_separated_polished_assemblies[org].append(assembly)

# Generate reference mapped scaffolds
ragtag_commands = ragtag_pipe.generate_ragtag_commands(org_separated_polished_assemblies, ragtag_results_directory,
                                                       threads, selected_reference_path_dictionary)

# Run the Ragtag commands with error handling
for command_all in ragtag_commands:
    try:
        subprocess.run(command_all[1], check=True)
    except subprocess.CalledProcessError as e:
        error_log.append(f'Ragtag error for specimen: {command_all[0]} with error code: {e.returncode}')

assemblies_final = defaultdict(list)
folders_ragtag = os.listdir(ragtag_results_directory)
directory_paths_ragtag = [(os.path.join(ragtag_results_directory, folder), folder)
                          for folder in folders_ragtag
                          if os.path.isdir(os.path.join(ragtag_results_directory, folder))]
for directory in directory_paths_ragtag:
    file_detected = os.path.isfile(os.path.join(directory[0], 'ragtag.scaffold.fasta'))
    if file_detected:
        id_sample = directory[1]
        polished_file = 'ragtag.scaffold.fasta'
        assemblies_final[organism_id_dictionary[id_sample]].append((id_sample, os.path.join(directory[0],
                                                                                            polished_file)))


# Create kaptive, resfinder, plasmidfinder, mlst, quast and prokka commands
kaptive_commands = kaptive_pipe.generate_kaptive_commands(assemblies_final,
                                                          kaptive_results_directory)
resfinder_commands = resfinder_pipe.generate_resfinder_commands(assemblies_final,
                                                                resfinder_results_directory)
plasmidfinder_commands = plasmid_finder_pipe.generate_plasmid_finder_commands(assemblies_final,
                                                                              plasmidfinder_results_directory)
mlst_commands = mlst_pipe.generate_mlst_commands(assemblies_final,
                                                 mlst_results_directory)
quast_commands = quast_pipe.generate_quast_commands(assemblies_final,
                                                    quast_results_directory, threads,
                                                    selected_reference_path_dictionary)
prokka_commands = prokka_pipe.generate_prokka_commands(assemblies_final,
                                                       prokka_results_directory)
virulence_commands = gen_find.generate_gen_db_commands(assemblies_final,
                                                       std_paths.virulence_db_path,
                                                       std_paths.virulence_db,
                                                       virulence_results_directory)

# Create commands array
commands = [quast_commands,
            resfinder_commands,
            plasmidfinder_commands,
            mlst_commands,
            prokka_commands,
            virulence_commands,
            kaptive_commands]

# Run the commands with error handling
for command_item in commands:
    for command in command_item:
        try:
            subprocess.run(command[1], check=True)
        except subprocess.CalledProcessError as e:
            error_log.append(f'{command[0]}: Error for command: {command[1]} with error code: {e.returncode}')

# Generate prokka results file for roary analysis
folders_prokka = os.listdir(prokka_results_directory)
directory_paths_prokka = [(os.path.join(prokka_results_directory, folder), folder)
                          for folder in folders_prokka
                          if os.path.isdir(os.path.join(prokka_results_directory, folder))]
for directory in directory_paths_prokka:
    folders_isolate = os.listdir(directory[0])
    org_results_roary = os.path.join(roary_results_directory, directory[1])
    if not os.path.exists(org_results_roary):
        os.makedirs(org_results_roary)
    isolate_paths = [(os.path.join(directory[0], folder), folder)
                     for folder in folders_isolate
                     if os.path.isdir(os.path.join(directory[0], folder))]
    for isolate in isolate_paths:
        file = os.path.join(isolate[0], f'{isolate[1]}.gff')
        if os.path.isfile(file):
            shutil.copy(file, org_results_roary)

# Generate roary commands
roary_commands = roary_pipe.generate_roary_commands(roary_results_directory, threads)


# Run roary commands
for command in roary_commands:
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        error_log.append(f'Roary command: {command} failed with error code: {e.returncode}')

# Generate dictionary for roary secondary analysis
# roary_dict = {}
# folders_roary = os.listdir(roary_results_directory)
# directory_paths_roary = [(os.path.join(roary_results_directory, folder), folder) for folder in folders_roary if os.path.isdir(os.path.join(roary_results_directory, folder))]
# for directory in directory_paths_roary:
#     folders_isolate = os.listdir(directory[0])
#     isolate_paths = [(os.path.join(directory[0], folder), directory[1]) for folder in folders_isolate if
#                      os.path.isdir(os.path.join(directory[0], folder))]
#     for isolate in isolate_paths:
#         gene_presence_csv = os.path.join(isolate[0], 'gene_presence_absence.csv')
#         alignment = os.path.join(isolate[0], 'core_gene_alignment.aln')
#         file_detected_csv = os.path.isfile(os.path.join(directory[0], 'gene_presence_absence.csv'))
#         if os.path.isfile(gene_presence_csv) and os.path.isfile(alignment):
#             roary_dict[directory[1]] = (gene_presence_csv, alignment)
#
# # Generate Fasttree and roary_plot commands
# fasttree_commands = fasttree_pipe.generate_fasttree_commands(roary_dict, fasttree_results_directory)
# roary_plot_commands = roary_plot_pipe.generate_roary_plot_commands(roary_dict, fasttree_results_directory)
#
# # Run Fasttree commands
# for command in fasttree_commands:
#     try:
#         subprocess.run(command[0], check=True, stdout=open(command[1], 'w'))
#     except subprocess.CalledProcessError as e:
#         error_log.append(f'Fasttree command: {command} failed with error code: {e.returncode}')
#
# # Run roary_plot commands
# for command in roary_plot_commands:
#     try:
#         subprocess.run(command[0], check=True)
#     except subprocess.CalledProcessError as e:
#         error_log.append(f'Roary_plot command: {command} failed with error code: {e.returncode}')
#

# Output data figures
data_handling.generate_figures_data(figure_results_directory,
                                    resfinder_results_directory,
                                    plasmidfinder_results_directory,
                                    mlst_results_directory,
                                    quast_results_directory,
                                    virulence_results_directory,
                                    kaptive_results_directory)

# Create parsnp folders for tree generation
parsnp_pipe.generate_evaluation_folder_parsnp(assemblies_final, parsnp_evaluation_directory)
folders_parsnp = os.listdir(parsnp_evaluation_directory)
directory_paths_parsnp = [(os.path.join(parsnp_evaluation_directory, folder), folder)
                          for folder in folders_parsnp
                          if os.path.isdir(os.path.join(parsnp_evaluation_directory, folder))]
for directory in directory_paths_parsnp:
    folders_isolate = os.listdir(directory[0])
    org_results_brig = os.path.join(brig_evaluation_directory, directory[1])
    if not os.path.exists(org_results_brig):
        os.makedirs(org_results_brig)
    isolate_paths = [(os.path.join(directory[0], folder), folder)
                     for folder in folders_isolate
                     if os.path.isfile(os.path.join(directory[0], folder))]
    for isolate in isolate_paths:
        file = isolate[0]
        if os.path.isfile(file):
            output_file = os.path.join(org_results_brig, f'{isolate[1]}')
            BRIG_file_creation.return_largest_sequence(file, output_file)

# Error log generation
for error in error_log:
    with open('test/error_log', 'a+') as f:
        f.seek(0, 2)
        f.write(f'{error}\n')








