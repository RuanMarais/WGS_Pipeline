from Bio import SeqIO


def return_largest_sequence(input_file, output_file):
    # Read and parse the multi-FASTA file
    sequences = list(SeqIO.parse(input_file, "fasta"))

    # Find the largest sequence
    largest_sequence = max(sequences, key=lambda seq: len(seq))

    # Write the largest sequence to a new FASTA file
    with open(output_file, "w") as out_handle:
        SeqIO.write(largest_sequence, out_handle, "fasta")

