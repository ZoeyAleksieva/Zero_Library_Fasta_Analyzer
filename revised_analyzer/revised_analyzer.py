import sample_qc
import sys

# The following code was lifted from the generator lesson of the Python for biologists course by ComputingLifeSci
# a generator function — looks like a normal function, but uses 'yield' instead of 'return'
def read_fasta_records(filename):
   with open(filename) as fasta_file:
       header = None
       sequence_lines = []
       for line in fasta_file:
           line = line.strip()
           if line.startswith(">"):
               if header is not None:
                   yield header, "".join(sequence_lines)   # produce one record, then pause here
               header = line[1:]
               sequence_lines = []
           else:
               sequence_lines.append(line)
       if header is not None:
           yield header, "".join(sequence_lines)   # yield the final record after the loop ends

# MAIN
if len(sys.argv) < 3:
    print("Usage: python main.py <fasta_file> qc_threshold")
    sys.exit(1)

#use try/catch to handle invalid input for the fasta file
try:
    fasta_file = sys.argv[1]
    with open(fasta_file, 'r') as f:
        pass  # Just to check if the file can be opened
except FileNotFoundError:
    print(f"Error: File '{fasta_file}' not found.")
    sys.exit(1)

#use try/catch to handle invalid input for the QC threshold
try:
    qc_threshold = float(sys.argv[2])
except ValueError:
    print("Error: QC threshold must be a number.")
    sys.exit(1)

fasta_file = sys.argv[1]
qc_threshold = float(sys.argv[2])  # Set the QC threshold for GC content
samples = []

for header, sequence in read_fasta_records(fasta_file):
   print(f"Sample ID: {header}, Sequence length: {len(sequence)} bases")
   curr = sample_qc.Sample(header=header, sequence=sequence)
   samples.append(curr)
   print(f"GC content: {curr.get_gc_content():.2f}")
   print(f"Reverse complement: {curr.get_reverse_complement()}")
   print(f"Longest orf seq: {curr.get_longest_ORF_seq()}")
   # Write the results to a file
   with open("results.txt", "a") as results_file:
       results_file.write(f"{header}\n")
       results_file.write(f"GC content: {curr.get_gc_content():.2f}\n")
       results_file.write(f"Reverse complement: {curr.get_reverse_complement()}\n")
       results_file.write(f"Longest orf seq: {curr.get_longest_ORF_seq()}\n\n")
       results_file.write("--------------------------------------------------\n\n")

#Use list comprehensions to filter the cohort down to samples passing a GC-content QC threshold
filtered_samples = [sample for sample in samples if sample.get_gc_content() >= qc_threshold]
for sample in filtered_samples:
    print(f"Sample {sample.sample_id} passed QC with GC content: {sample.get_gc_content():.2f}")