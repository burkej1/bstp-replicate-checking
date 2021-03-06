from __future__ import print_function
from cyvcf2 import VCF
import sys
import os
import re
import numpy as np
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Takes a folder location containing vcfs with paired " \
                                                 "samples, calculates concordance and prints to a given " \
                                                 "output file. Assumes vcfs to be named .merged.BSID.vcf")
    parser.add_argument("-i", "--inputfolder", help="Folder containing vcf files.", required=True)
    parser.add_argument("-o", "--outputtsv", help="Location to write the summary tsv.", required=True)
    return vars(parser.parse_args())


#vcf_folder_path = "/scratch/vh83/projects/brastrap/run1-7_replicate_pipeline/variants/gatk"

class Concordance(object):
    """Class to keep track of concordance metrics between two samples."""
    def __init__(self, samples):
        # Meanings of gt_types values
        genotype_dictionary = {0: "hom_ref", 1: "het", 2: "unknown", 3: "hom_alt"}
        # Keeping track of concordance
        self.concordant_calls = 0.0
        self.discordant_calls = 0.0
        # Sets for discordance reason determination
        self.het_hom = set([1, 3])
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
        else:
            print("No discordance reason found")
            print(set_genotypes)

    def calculate_concordance(self):
        try:
            concordance = self.concordant_calls / (self.concordant_calls + self.discordant_calls)
        except ZeroDivisionError:
            concordance = 0
        return concordance


def get_individual_vcfs(folder_path):
    folder_contents = os.listdir(folder_path)
    vcf_list = [folder_path + '/' + vcf for vcf in folder_contents if re.match(".+\.merged\.BS\d\d\d\d\d\d\.vcf$", vcf)]
    return vcf_list


def snp_filter(variant):
    """SNP filtering thresholds, only prints a warning if a field that shouldn't be missing is missing. MQRankSum
    and ReadPosRankSum seem to not be calculated relatively often."""
    fail = False
    missing_metric = False
    if variant.QUAL < 30.0:
        fail = True
    if variant.INFO.get('DP') < 10.0 and variant.INFO.get('DP') is not None:
        fail = True
    elif variant.INFO.get('DP') is None:
        print("DP missing")
        missing_metric = True
    if variant.INFO.get('QD') < 2.0 and variant.INFO.get('QD') is not None:
        fail = True
    elif variant.INFO.get('QD') is None:
        print("QD missing")
        missing_metric = True
    if variant.INFO.get('MQ') < 30.0 and variant.INFO.get('MQ') is not None:
        fail = True
    elif variant.INFO.get('MQ') is None:
        print("MQ missing")
        missing_metric = True
    if variant.INFO.get('SOR') > 3.0 and variant.INFO.get('SOR') is not None:
        fail = True
    elif variant.INFO.get('SOR') is None:
        print("SOR missing")
        missing_metric = True
    if variant.INFO.get('MQRankSum') < -12.5 and variant.INFO.get('MQRankSum') is not None:
        fail = True
    elif variant.INFO.get('MQRankSum') is None:
        missing_metric = True
    if variant.INFO.get('ReadPosRankSum') < -8.0 and variant.INFO.get('ReadPosRankSum') is not None:
        fail = True
    elif variant.INFO.get('ReadPosRankSum') is None:
        missing_metric = True
    return fail


def indel_filter(variant):
    """Indels filtering thresholds, only prints a warning if a field that shouldn't be missing is missing. 
    ReadPosRankSum seems to not be calculated relatively often."""
    fail = False
    missing_metric = False
    if variant.QUAL < 30.0:
        fail = True
    if variant.INFO.get('DP') < 10.0 and variant.INFO.get('DP') is not None:
        fail = True
    elif variant.INFO.get('DP') is None:
        print("DP missing")
        missing_metric = True
    if variant.INFO.get('QD') < 2.0 and variant.INFO.get('QD') is not None:
        fail = True
    elif variant.INFO.get('QD') is None:
        print("QD missing")
        missing_metric = True
    if variant.INFO.get('ReadPosRankSum') < -20.0 and variant.INFO.get('ReadPosRankSum') is not None:
        fail = True
    elif variant.INFO.get('ReadPosRankSum') is None:
        missing_metric = True
    return fail


def calculate_concordance(vcf_path, outputfile):
    """Calculates concordance for vcf samples in an individual vcf."""
    BSID = re.search(".+\.merged\.(BS\d\d\d\d\d\d)\.vcf", vcf_path).group(1)
    concordance_metrics = Concordance(BSID)
    vcf = VCF(vcf_path)

    # Counting the variants filtered
    filtered_snps = 0
    filtered_indels = 0

    # Keeping track of the total variants and plate number on a per-sample basis
    genotype_dictionary = {0: "hom_ref", 1: "het", 2: "unknown", 3: "hom_alt"}
    total_vars = {}
    samples = vcf.samples
    for sample in samples:
        platenumber = re.search("_(\d\d\d\d)-", sample).group(1)
        total_vars[sample] = {"hom_ref": 0, "hom_alt": 0, "unknown": 0, "het": 0, "plate": platenumber}

    if total_vars[samples[0]]["plate"] == total_vars[samples[1]]["plate"]:
        inter_intra = "Intraplate"
    else:
        inter_intra = "Interplate"

    # Iterating through the variants in the vcf and calculating sample concordance
    for variant in vcf:
        # Statistics pre-filtering (looking for strange amounts of noise)
        for n, gt in enumerate(variant.gt_types):
            sample = samples[n]
            genotype = genotype_dictionary[gt]
            total_vars[sample][genotype] += 1
        # Filtering, skip if the variant fails filters
        if variant.is_snp:
            if variant.FILTER is not None:
                filtered_snps += 1
                continue
        elif variant.is_indel:
            if variant.FILTER is not None:
                filtered_indels += 1
                continue
        else:
            print("Not a SNP or Indel, might want to look into this variant.")
            print(variant)
        # Checking concordance (post-filtering)
        genotypes = variant.gt_types
        set_genotypes = set(genotypes)
        if genotypes[0] == genotypes[1]:
            concordance_metrics.concordant_calls += 1.0
        else:
            concordance_metrics.discordant_calls += 1.0
            concordance_metrics.check_discordance_reason(set_genotypes)
    outputline = "{bsid}\t{samples}\t{interorintra}\t{concordance}%\t{filteredindels}\t{filteredsnps}\t" \
                 "{concordantcalls}\t{discordantcalls}\t{het_hom}\t{hom_var}\t{call_nocall}\t" \
                 "{unf_het}\t{unf_homref}\t{unf_homalt}\t{unf_nocall}\n" \
                 .format(bsid=BSID, 
                         samples="{s1}, {s2}".format(s1=samples[0], s2=samples[1]),
                         interorintra=inter_intra, 
                         concordance=str(round(concordance_metrics.calculate_concordance() * 100, 2)), 
                         filteredindels=filtered_indels, 
                         filteredsnps=filtered_snps, 
                         concordantcalls=concordance_metrics.concordant_calls,
                         discordantcalls=concordance_metrics.discordant_calls, 
                         het_hom=concordance_metrics.discordance_reasons["het_hom"], 
                         hom_var=concordance_metrics.discordance_reasons["ref_var"], 
                         call_nocall=concordance_metrics.discordance_reasons["no_call"], 
                         unf_het="{s1}//{s2}".format(s1=total_vars[samples[0]]["het"], 
                                                    s2=total_vars[samples[1]]["het"]), 
                         unf_homref="{s1}//{s2}".format(s1=total_vars[samples[0]]["hom_ref"], 
                                                       s2=total_vars[samples[1]]["hom_ref"]), 
                         unf_homalt="{s1}//{s2}".format(s1=total_vars[samples[0]]["hom_alt"], 
                                                       s2=total_vars[samples[1]]["hom_alt"]), 
                         unf_nocall="{s1}//{s2}".format(s1=total_vars[samples[0]]["unknown"], 
                                                       s2=total_vars[samples[1]]["unknown"]))
    outputfile.write(outputline)


def main():
    args = parse_args()
    if args["inputfolder"].endswith('/'):
        args["inputfolder"] = args["inputfolder"].rstrip('/')
    vcf_list = get_individual_vcfs(args["inputfolder"])
    header = "BSID\tSamples\tInter- or Intra-\tConcordance\tFiltered Indels (number removed)\t" \
             "Filtered SNPs (number removed)\tConcordant Calls\tDiscordant Calls\t" \
             "Heterozygous vs. Homozygous Alt\tHomozygous Reference vs Variant\tCall vs. No Call\t" \
             "Unfiltered Het Calls\tUnfiltered Hom Ref Calls\t Unfiltered Hom Alt Calls\t" \
             "Unfiltered No Calls\n"
    with open(args["outputtsv"], 'w') as outputfile:
        outputfile.write(header)
        for vcf in vcf_list:
            calculate_concordance(vcf, outputfile)


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
