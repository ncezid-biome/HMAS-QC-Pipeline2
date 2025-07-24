#!/usr/bin/env python

import yaml
import sys
import re
from collections import defaultdict
import os

def parse_header_info(header_info_str):
    """Parses the multiqc_header string into a dictionary."""
    # Remove brackets and split key-value pairs
    header_info_str = header_info_str.strip("[]")
    key_value_pairs = re.findall(r'([^:,]+):([^,]+)', header_info_str)

    # Convert to dictionary, stripping extra spaces
    return {k.strip(): v.strip() for k, v in key_value_pairs}

def parse_software_version(software_versions):
    # Load the original versions.yml
    with open(software_versions) as f:
        raw_versions = yaml.safe_load(f)

    # Flatten the structure and accumulate versions in each section
    software_versions = defaultdict(dict)
    for section, versions in raw_versions.items():
        for k, v in versions.items():
            # If the section already has this key, append it, otherwise add it
            if k in software_versions[section]:
                software_versions[section][k] = v
            else:
                software_versions[section][k] = v

    # Flatten to a dictionary with all versions
    flat_versions = {}
    for section, versions in software_versions.items():
        flat_versions.update(versions)

    return flat_versions


def update_multiqc_config(yaml_file, header_info_str, software_versions, primer_file):
    """Appends parsed header info and software versions to multiqc_config.yaml while maintaining structure."""
    
    # Parse header_info
    header_info_dict = parse_header_info(header_info_str)
    software_version_dict = parse_software_version(software_versions)

    # Load existing YAML file
    with open(yaml_file, "r") as f:
        config = yaml.safe_load(f) or {}

    # create a new 'report_header_info' even if it exists already
    config["report_header_info"] = []

    # Append new header info without duplication
    for key, value in header_info_dict.items():
        entry = {key: str(value)}
        config["report_header_info"].append(entry)
    
    # add the primer / oligo information
    config["report_header_info"].append({"oligo file": os.path.basename(primer_file)})

    # Add the software_versions section (making sure to not have dashes)
    config["software_versions"] = software_version_dict

    # Write back to YAML file
    with open(yaml_file, "w") as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)


if __name__ == "__main__":
    yaml_file = sys.argv[1]  # First argument: path to multiqc_config.yaml
    header_info_str = sys.argv[2]  # Second argument: multiqc_header as a string
    software_versions = sys.argv[3] # 3rd argument: software_version.yaml 
    primer_file = sys.argv[4] # 4th argument: primer file (oligo file)  
    update_multiqc_config(yaml_file, header_info_str, software_versions, primer_file)
