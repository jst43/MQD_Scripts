"""Moves to target directory, and for all fastq.gz files a Sample instance is created.
These objects are used to ensure correct input names for files, as well as correct
output names. Other QC checks ensure that each biological sample has 4 files (2 index,
2 read) and defines which fastq type ('R1','R2','R3','I1','I2') is associated with index
for deduplication, and which types correspond to Read 1 and Read 2."""

# PACKAGES
import os
import re
import argparse
import gzip
from Bio import SeqIO

# CLASSES and FUNCTIONS
class Sample(object):
  """Object to hold information from fastq files,
  allowing grouping by biological sample and QC checking."""
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
  def parse_sample_name(sample_name):
    dup_match = re.search('dup', sample_name.lower())
    if dup_match is not None:
      dup_span = dup_match.span()
      sample_name = sample_name[0:dup_span[0] - 1] +\
                        "_Dup" + sample_name[dup_span[1]:]
      try:
        int(sample_name[dup_span[1]:])
      except ValueError:
        print('Unexpected value for dup number')
    return sample_name

  @staticmethod
  def parse_sample_filestring(string):
    """Assumes the input filestring is of the form
    [SampleName](-Dup[N])_[Sample_Num]_[Lane]_[fastq_type]_001[file_extension]
    Where SampleNum is of the form 'S[1-9]+'
    and Lane is of the form 'L00[1-9]'"""
    sample_name, suffix = string.split("_S")
    sample_name = Sample.parse_sample_name(sample_name)
    sample_list = suffix.split("_")
    try:
      assert len(sample_list) == 4
    except AssertionError:
      raise Exception('Filestring is in unexpected format.')
    # Check that file extension is 'fastq.gz'
    file_extension = sample_list[-1].replace("001", "")
    try:
      assert file_extension == ".fastq.gz"
    except AssertionError:
      raise Exception('File is not a gzipped fastq.')
    sample_num = "S" + sample_list[0]
    _, lane, fastq_type, _ = sample_list
    return sample_name, sample_num, lane, fastq_type, file_extension

  def filename_output(self):
    return self.name + "_" + self.lane + "_" + self.type + self.extension

  def get_read_length(self):
    file_open = gzip.open(self.inputfile, 'rt')
    record = next(SeqIO.parse(file_open, 'fastq'))
    read_length = len(record)
    file_open.close()
    return read_length

def read_samples():
  fastq_samples = list()
  for file in os.listdir(os.getcwd()):
    if file.endswith('.fastq.gz'):
      file_sample = Sample.from_filename(file)
      fastq_samples.append(file_sample)
  return fastq_samples

def check_samples(fastq_samples):
  """For list of Samples, count number of files for each Biological Sample
  and check read length. Index files have outputfile written to list,
  But for read_files (the 151bp fastqs) there are 2 sets of files, 'R1' and 'R2'.
  Because of demultiplexing, 'R1' and 'R2' files may not be named R1 and R2 in Sample.type
  So read files have entire Sample object appended to list for further sorting."""
  index_read_length = 10
  std_read_length = 151
  sample_names_count = dict()
  index_type = set()
  read_type_unique = set()
  for sample in fastq_samples:
    if sample.name in sample_names_count:
      sample_names_count[sample.name] += 1
    else:
      sample_names_count[sample.name] = 1
    #Sort by read_length"
    if sample.read_length == index_read_length:
      index_type.add(sample.type)
    elif sample.read_length == std_read_length:
      read_type_unique.add(sample.type)
  #Check that each sample has 4 files"
  if all(value == 4 for value in sample_names_count.values()):
    print("All samples have 4 fastq.gz files")
  else:
    print(sample_names_count.keys())
    print(sample_names_count.values())
    raise Exception('Not all samples have 4 files')
  return index_type, read_type_unique

def sort_samples(fastq_samples, index_type, read_type_unique):
  """Sort read_files files by name"""
  try:
    assert len(read_type_unique) == 2
  except AssertionError:
    raise Exception('There should only be two classes of read fastq, %i found' % len(read_type_unique))
  try:
    assert len(index_type) == 1
  except AssertionError:
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
    handle.write('Index type:{}\n'.format(list(index_type)[0]))
    handle.write('R1 type:{}\n'.format(read_type_sorted[0]))
    handle.write('R2 type:{}\n'.format(read_type_sorted[1]))

def main():
    fastq_samples = read_samples()
    index_type, read_type_unique = check_samples(fastq_samples)
    sort_samples(fastq_samples, index_type, read_type_unique)

# FLAG ARGUMENTS
parser = argparse.ArgumentParser(description='Checks for filestring format and identifies index files')
parser.add_argument('-f', help='Specifies directory', type=str)
args = parser.parse_args()
filepath = args.f
os.chdir(filepath)

#SCRIPT
if __name__ == '__main__':
    main()