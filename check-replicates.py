from __future__ import print_function
from cyvcf2 import VCF
import sys
import os
import re
import numpy as np

#input_vcf_path = sys.argv[1]

vcf_folder_path = "/scratch/vh83/projects/brastrap/run1-7_replicate_pipeline/variants/gatk"
if vcf_folder_path.endswith('/'):
    vcf_folder_path = vcf_folder_path.rstrip('/')


class Concordance(object):
    """Class to keep track of concordance metrics between two samples."""
    def __init__(self, samples):
        # Meanings of gt_types values
        genotype_dictionary = {0: "hom_ref", 1: "het", 2: "unknown", 3: "hom_alt"}
        # Keeping track of concordance
        self.concordant_calls = 0.0
        self.discordant_calls = 0.0
        # Sets for discordance reason determination
        self.het_hom = set([2, 3])
        self.ref_var_het = set([0, 1])
        self.ref_var_hom = set([0, 3])
        self.discordance_reasons = {"no_call": 0, "het_hom": 0, "ref_var": 0}

    def check_discordance_reason(self, set_genotypes):
        # Checks reason for discordance and updates counts
        if set_genotypes == self.het_hom:
            self.discordance_reasons["het_hom"] += 1
        elif set_genotypes == self.ref_var_hom or set_genotypes == self.ref_var_het:
            self.discordance_reasons["ref_var"] += 1
        elif 2 in set_genotypes:
            self.discordance_reasons["no_call"] += 1

    def calculate_concordance(self):
        try:
            concordance = self.concordant_calls / (self.concordant_calls + self.discordant_calls)
        except ZeroDivisionError:
            concordance = 0
        return concordance


def get_individual_vcfs(folder_path):
    folder_contents = os.listdir(folder_path)
    vcf_list = [folder_path + '/' + vcf for vcf in folder_contents if re.match(".+\.BS\d\d\d\d\d\d\.vcf$", vcf)]
    return vcf_list


def snp_filter(variant):
    fail = False
    qual_f = 30.0
    qd_f = 2.0
    dp_f = 10.0
    mq_f = 30.0
    sor_f = 3.0
    mqrs_f = -12.5
    rprs_f = -8.0

    return fail


def indel_filter(variant):
    fail = False
    return fail


def calculate_concordance(vcf_path):
    """Calculates concordance for vcf samples in an individual vcf."""
    BSID = re.search(".+\.(BS\d\d\d\d\d\d)\.vcf", vcf_path).group(1)
    concordance_metrics = Concordance(BSID)
    vcf = VCF(vcf_path)
    for variant in vcf:
        # Filtering, skip if the variant fails filters
        if is_snp:
            if snp_filter(variant) is True:
                continue
        if is_indel:
            if indel_filter(variant) is True:
                continue
        # Checking concordance
        genotypes = variant.gt_types
        set_genotypes = set(genotypes)
        if genotypes[0] == genotypes[1]:
            concordance_metrics.concordant_calls += 1.0
        else:
            concordance_metrics.discordant_calls += 1.0
            concordance_metrics.check_discordance_reason(set_genotypes)
    print(concordance_metrics.calculate_concordance())
    print("Concordant calls: {conc}\n" \
          "Discordant calls: {disc}".format(conc=concordance_metrics.concordant_calls, 
                                            disc=concordance_metrics.discordant_calls))
    print(concordance_metrics.discordance_reasons)


def main():
    vcf_list = get_individual_vcfs(vcf_folder_path)
    for vcf in vcf_list:
        calculate_concordance(vcf)


if __name__ == "__main__":
    main()











#def calculate_concordance(duplicate_dictionary, vcf_path):
#    """For each duplicate pair, iterates through all variants in the vcf and calculates concordance/discordance
#    along with reasons for discordance."""
#    for ID in duplicate_dictionary:
#        concordance_metrics = Concordance(duplicate_dictionary[ID])
#        vcf = VCF(vcf_path)
#        # Subsets the vcf to the samples of interest
#        vcf.set_samples(duplicate_dictionary[ID])
#        for variant in vcf:
#            # Skip if variant has been filtered
#            if variant.FILTER != None:
#                continue
#            # Skip if both samples are hom_ref at this location
#            if variant.num_hom_ref == 2 or variant.num_unknown == 2:
#                continue
#            genotypes = variant.gt_types
#            set_genotypes = set(genotypes)
#            if genotypes[0] == genotypes[1]:
#                concordance_metrics.concordant_calls += 1.0
#            else:
#                concordance_metrics.discordant_calls += 1.0
#                concordance_metrics.check_discordance_reason(set_genotypes)
#        print(concordance_metrics.calculate_concordance())
#        print("Concordant calls: {conc}\n" \
#              "Discordant calls: {disc}".format(conc=concordance_metrics.concordant_calls, 
#                                                disc=concordance_metrics.discordant_calls))
#        print(concordance_metrics.discordance_reasons)

#def retrieve_duplicate_ids(vcf_path):
#    # Dictionary to store full sample IDs with their BS ID as a key
#    sample_dictionary = {}
#    for sample in VCF(vcf_path).samples:
#        if not re.match(".+BS\d\d\d\d\d\d.+", sample):
#            print("No BS ID found for sample: {sample}".format(sample=sample))
#            continue
#        sampleid = re.search("BS\d\d\d\d\d\d", sample).group(0)
#        if sampleid in sample_dictionary:
#            sample_dictionary[sampleid].append(sample)
#        else:
#            sample_dictionary[sampleid] = [sample]
#    # Adding duplicates only to the duplicate dictionary
#    duplicate_dictionary = {}
#    for key in sample_dictionary:
#        if len(sample_dictionary[key]) > 1:
#            duplicate_dictionary[key] = sample_dictionary[key]
#        else:
#            print("Duplicate not found for BSID: {}".format(key))
#    return duplicate_dictionary
