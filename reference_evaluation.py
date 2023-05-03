import os
import re


def extract_avg_identity(file_path):
    with open(file_path, 'r') as file:
        content = file.read()
        avg_identity_pattern = r"AvgIdentity\s+\d+\.\d+\s+(\d+\.\d+)"
        one_to_one_section_pattern = r"\[Alignments\][\s\S]*?1-to-1[\s\S]*?AvgIdentity\s+\d+\.\d+\s+(\d+\.\d+)"

        one_to_one_section = re.search(one_to_one_section_pattern, content)
        if one_to_one_section:
            avg_identity_match = re.search(avg_identity_pattern, one_to_one_section.group(0))
            if avg_identity_match:
                return float(avg_identity_match.group(1))
    return None


def evaluate_average_identity(directory_path):
    directory = directory_path
    report_files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('mummer.report')]

    identities = []
    for report_file in report_files:
        avg_identity = extract_avg_identity(report_file)
        if avg_identity is not None:
            identities.append(avg_identity)

    if identities:
        overall_avg_identity = sum(identities) / len(identities)
        return overall_avg_identity
    else:
        return None

