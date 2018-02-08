# Code for checking sample replicates

check-replicates.py takes a folder location as input and looks for individually genotyped replicate pairs in the following format SOMESAMPLENAME.BS123456.vcf
It then calculates concordance and a number of other metrics and writes it all to an outputfile in the current working directory. Currently automatically named to replicate_metrics.tsv

generate_filelist.py takes a list of fastqs as input, extracts and pairs replicates and writes to an output file. This can then be used to run the replicate genotyping pipeline.
