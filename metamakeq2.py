## Qiime2 metadata maker
# Made by Group 5 for BIOT 670I
# This script will read in a MIMARKS-style metadata file and select the user-specified columns only 
# for use in Qiime2 as a metadata file. Mostly intended for use with data from SRA. 
import csv


print("**For maximum compatibility the names of fields should be no more than 10 characters long.")
print("**If you encounter issues with Qiime an explicit type declaration line may be needed")

filename = input("Please enter the name of the metadata file, with extension: ")
metafields = input("Please enter a space-separated list of the exact column names you would like to include: ")
newfields = input("Please enter a space-separated list of new names for the columns for use in Qiime2 \n(Note: #SampleID is already included, skip it): ")

# extract new Qiime-compatible metadata file column names from user's list
headerlist = ["#SampleID"]
columns = newfields.split()
for i in columns:
    headerlist.append(i)
    
headerline = "\t".join(headerlist)
print(f"{headerline}\n")

# convert metafields into an f string of the right size.
metalist = metafields.split()

#
outlist = []

# block to open and read file
with open(filename, "r") as metasra, open("metadataq2.tsv", "w") as q2meta:
    reader = csv.DictReader(metasra, delimiter=',')
    #next(reader)
    q2meta.write(f"{headerline}\n")
    #q2meta.write('#SampleID\tgroup\tsex\n') # example header line
    #q2meta.write('#q2:types\tcategorical\tcategorical\n') # example types line # (uncomment and modify if needed)
    for line in reader:
        #q2meta.write(f"{line['host_subject_id']}\t{line['Host_disease']}\t{line['host_sex']}\t{line['host_body_mass_index']}\n") # columns used for dataset.
        for i in metalist:
            q2meta.write(f"{line[i]}\t")
        q2meta.write("\n")
        


        

