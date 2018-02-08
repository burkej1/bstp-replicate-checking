import sys
import re

input_file_path = sys.argv[1]
output_file_path = "pipeline_formatted_filelist.txt"

id_dictionary = {}

with open(input_file_path, 'r') as input_file, open(output_file_path, 'w') as output_file:
    # Reading filenames into a dictionary
    for line in input_file:
        stripline = line.rstrip('\n')
        if "_R2_" in stripline:  # Save time by only looking at R1 fastqs
            continue

        # Check for a match then add filenames to the dictionary
        bs_id_result = re.search("-(BS\d\d\d\d\d\d).+R1", stripline)
        if bs_id_result:  # If there's a BS ID add the filename to the dictionary under the ID
            bs_id = bs_id_result.group(1)
            if bs_id in id_dictionary:
                id_dictionary[bs_id].append(stripline)
            else:
                id_dictionary[bs_id] = [stripline]

    replicate_dictionary = {}
    for ID in id_dictionary:
        if len(id_dictionary[ID]) > 1:
            is_real_replicate = False
            for filename in id_dictionary[ID]:
                # Making sure at least one of the files is tagged as a replicate
                if "_R_" in filename:
                    is_real_replicate = True
            if is_real_replicate:
                replicate_dictionary[ID] = id_dictionary[ID]

    # replicate_numbers = [len(replicate_dictionary[bsid]) for bsid in replicate_dictionary]
    # for ID in replicate_dictionary:
    #     print ID
    #     for filename in replicate_dictionary[ID]:
    #         if filename in fails:
    #             print filename + "\t\tFAIL"
    #         else:
    #             print filename
    #     print '\n'

    # Writing replicates to an output file to cat to the pipeline, or create symlinks
    for bsid in replicate_dictionary:
        for replicate in replicate_dictionary[bsid]:
            output_file.write(replicate + '\n')
            output_file.write(re.sub("R1", "R2", replicate) + '\n')



# # Controls, removed for now
#         if re.match(".*run7.+ABCFS.+H12.+", stripline):  # Checking for positive controls with BS IDs
#             id_dictionary["positive_controls"].append(stripline)
#             continue
# 
#         else:  # If there's no BS ID check if it's a control
#             pc_result = re.search(".+X4336.+", stripline)
#             if pc_result:
#                 id_dictionary["positive_controls"].append(stripline) 
# "positive_controls": []
#     for control in id_dictionary["positive_controls"]:
#         if control in fails:
#             print "FAIL"
#             print control
#         output_file.write(control + '\n')
#         output_file.write(re.sub("R1", "R2", control) + '\n')
