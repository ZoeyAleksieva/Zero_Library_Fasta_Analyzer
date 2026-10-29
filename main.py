import sys
import analyze_fasta

if len(sys.argv) < 2:
    print("Usage: python main.py <fasta_file>")
    sys.exit(1)

#TODO: Check if it's a valid fasta file

fasta_file = sys.argv[1]
analyze_fasta(fasta_file)

