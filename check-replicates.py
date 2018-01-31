from cyvcf2 import VCF
import sys
import re
import random  # For testing (no replicates in the testing vcf so going to randomly pair sample IDs)
import numpy as np

input_vcf_path = sys.argv[1]

def retrieve_duplicate_ids(vcf_path):
    count = 0  # TESTING
    # Dictionary to store full sample IDs with their BS ID as a key
    sample_dictionary = {}
    for sample in VCF(vcf_path).samples:
        if not re.match(".+BS\d\d\d\d\d\d.+", sample):
            print "No BS ID found for sample: {sample}".format(sample=sample)
            continue
        sampleid = re.search("BS\d\d\d\d\d\d", sample).group(0)
        if sampleid in sample_dictionary:
            continue  # TESTING
            sample_dictionary[sampleid].append(sample)
        else:
            sample_dictionary[sampleid] = [sample]
            sample_dictionary[sampleid].append(random.choice(VCF(vcf_path).samples))  # TESTING
            count += 1  # TESTING
            if count == 20:  # TESTING
                break  # TESTING

    # Adding duplicates only to the duplicate dictionary
    duplicate_dictionary = {}
    for key in sample_dictionary:
        if len(sample_dictionary[key]) > 1:
            duplicate_dictionary[key] = sample_dictionary[key]

    return duplicate_dictionary


def calculate_concordance(duplicate_dictionary, vcf_path):
    print "Starting iteration"
    for ID in duplicate_dictionary:
        vcf = VCF(vcf_path)
        # Subsets the vcf to the samples of interest
        vcf.set_samples(duplicate_dictionary[ID])
        for variant in vcf:
            # Selects variant sites
            if not variant.num_hom_ref == 2:
                pass
    print "All complete (20)."


def main():
    duplicate_ids = retrieve_duplicate_ids(input_vcf_path)
    calculate_concordance(duplicate_ids, input_vcf_path)


if __name__ == "__main__":
    main()











