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
    """For each duplicate pair, iterates through all variants in the vcf and calculates concordance/discordance
    along with reasons for discordance."""
    # Meanings of gt_types values
    genotype_dictionary = {0: "hom_ref", 1: "het", 2: "unknown", 3: "hom_alt"}
    # Sets for discordance reason determination
    het_hom = set([2, 3])
    ref_var_het = set([0, 1])
    ref_var_hom = set([0, 3])
    for ID in duplicate_dictionary:
        concordant_calls = 0.0
        discordant_calls = 0.0
        discordance_reasons = {"no_call": 0, "het_hom": 0, "ref_var": 0}
        vcf = VCF(vcf_path)
        # Subsets the vcf to the samples of interest
        vcf.set_samples(duplicate_dictionary[ID])
        for variant in vcf:
            # Skip if both samples are hom_ref at this location
            if variant.num_hom_ref == 2:
                continue
            genotypes = variant.gt_types
            set_genotypes = set(genotypes)
            if genotypes[0] == genotypes[1]:
                concordant_calls += 1.0
            else:
                discordant_calls += 1.0
                if 2 in set_genotypes:
                    discordance_reasons["no_call"] += 1
                elif set_genotypes == het_hom:
                    discordance_reasons["het_hom"] += 1
                elif set_genotypes == ref_var_het or set_genotypes == ref_var_hom:
                    discordance_reasons["ref_var"] += 1
                
        concordance = concordant_calls / (concordant_calls + discordant_calls)
        print concordance
        print discordance_reasons


def main():
    duplicate_ids = retrieve_duplicate_ids(input_vcf_path)
    calculate_concordance(duplicate_ids, input_vcf_path)


if __name__ == "__main__":
    main()











