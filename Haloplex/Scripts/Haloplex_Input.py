import os
import argparse
import gzip
from Bio import SeqIO

# CONSTANTS
index_read_length = 10
std_read_length = 151

# CLASSES and FUNCTIONS
class Sample(object):
    def __init__(self, filename, sample_name, sample_num, lane, fastq_type, extension):
        self.inputfile = filename
        self.name = sample_name
        self.num = sample_num
        self.lane = lane
        self.type = fastq_type
        self.extension = extension
        self.outputfile = self.filename_output()
        self.read_length = self.get_read_length()

    @classmethod
    def from_filename(cls, filename):
        sample_name, sample_num, lane, fastq_type, extension = cls.parse_sample_filestring(filename)
        return cls(filename, sample_name, sample_num, lane, fastq_type, extension)

    @staticmethod
    def parse_sample_filestring(str):
        sample_list = str.split("_")
        try:
            assert len(sample_list) == 5 or len(sample_list) == 6
        except:
            raise Exception('Filestring is in unexpected format.')
        file_extension = sample_list[-1].replace("001", "")
        try:
            assert file_extension == ".fastq.gz"
        except:
            raise Exception('File is not a gzipped fastq.')
        if len(sample_list) == 6:
            duplicate_name = sample_list[1].lower()
            try:
                assert duplicate_name[0:3] == "dup"
                assert isinstance(int(duplicate_name[-1]), int)
            except:
                raise Exception("Unexpected format for 'dup'.")
            sample_name = sample_list[0] + "-" + duplicate_name
            _, _, sample_num, lane, fastq_type, _ = sample_list
        else:
            sample_name = sample_list[0]
            _, sample_num, lane, fastq_type, _ = sample_list
        return (sample_name, sample_num, lane, fastq_type, file_extension)

    def filename_output(self):
        return self.name + "_" + self.lane + "_" + self.type + self.extension

    def get_read_length(self):
        file_open = gzip.open(self.inputfile, 'rt')
        record = next(SeqIO.parse(file_open, 'fastq'))
        read_length = len(record)
        file_open.close()
        return read_length

# FLAG ARGUMENTS
parser = argparse.ArgumentParser(description='Checks for filestring format and identifies index files')
parser.add_argument('f', help='Specifies directory', type=str)
args = parser.parse_args()
filepath = args.f
os.chdir(filepath)

# SCRIPT
# Read in fastq files as Samples
fastq_samples = list()
for file in os.listdir(os.getcwd()):
    if file.endswith('.fastq.gz'):
        file_sample = Sample.from_filename(file)
        fastq_samples.append(file_sample)

# For list of Samples, count number of files for each Biological Sample
# and check read length. Index files have outputfile written to list,
# but for read_files (the 151bp fastqs) there are 2 sets of files, 'R1' and 'R2'.
# Because of demultiplexing, 'R1' and 'R2' files may not be named R1 and R2 in Sample.type
# So read files have entire Sample object appended to list for further sorting.
sample_names_count = dict()
index_type = set()
read_type_unique = set()
for sample in fastq_samples:
    if sample.name in sample_names_count:
        sample_names_count[sample.name] += 1
    else:
        sample_names_count[sample.name] = 1
    #Sort by read_length
    if sample.read_length == index_read_length:
        index_type.add(sample.type)
    elif sample.read_length == std_read_length:
        read_type_unique.add(sample.type)


if all(value == 4 for value in sample_names_count.values()):
    print("All samples have 4 fastq.gz files")
else:
    print(sample_names_count.keys())
    print(sample_names_count.values())
    raise Exception('Not all samples have 4 files')

#Sort read_files files by name
try:
    assert len(read_type_unique) == 2
except:
    raise Exception('There should only be two classes of read fastq, %i found' % len(read_type_unique))

try:
    assert len(index_type) == 1
except:
    raise Exception('There should only be one class of read fastq, %i found' % len(index_type))
#Find which read_type has the smallest number, that read type will be 'Read 1'
# the other will be 'Read 2'
read_type_sorted = sorted(list(read_type_unique))

inputnames = [i.inputfile+'\n' for i in fastq_samples]
outputnames = [i.outputfile+'\n' for i in fastq_samples]

with open("original_names.txt", 'w') as handle:
    handle.writelines(sorted(inputnames))

with open("replacement_names.txt", 'w') as handle:
    handle.writelines(sorted(outputnames))

with open('fastq_types.txt', 'w') as handle:
    handle.write('Index type:{}\n'.format(index_type))
    handle.write('R1 type:{}\n'.format(read_type_sorted[0]))
    handle.write('R2 type:{}\n'.format(read_type_sorted[1]))