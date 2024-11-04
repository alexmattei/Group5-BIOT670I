### manifestmake q2
### Written by group 5 for BIOT 670I

import csv

print("**For maximum compatibility the names of fields should be no more than 10 characters long.")
print("**Use Bash pwd command to get the absolute filepath to the trimmed fastq seqs, and copy it. \n** Do not forget the forward slash (/) at the end when typing in path")

filename = input("Please enter the name of the Run Selector metadata file to parse, with extension: ")
sampleid = input("Please enter the exact column name to be used as a sample id: ")
dirpath = input("Please enter the absolute file path to the import directory with trimmed fastq files: ")

# block to open and read file
with open(filename, "r") as metasra, open("manifest-q2.tsv", "w") as q2manifest:
    reader = csv.DictReader(metasra, delimiter=',')
    #next(reader)
    q2manifest.write('sample-id\tabsolute-filepath\n')
    for line in reader:
        #q2manifest.write(f"{line['biospecimen_repository_sample_id']}\t{line['Run']}.fastq\n") # example of a common column name.
        q2manifest.write(f"{line[sampleid]}\t{dirpath}{line['Run']}.fastq\n")


