
# Zero-Library FASTA Analyzer

This repository contains two versions of a "zero-library" FASTA analyzer, built using concepts from introductory Python.

## 1. Simple Analyzer

The simple analyzer (`basic_analyzer.py` and `analyze_fasta.py`) uses fundamental Python concepts, including strings, dictionaries, loops, functions, and file handling.

- **`basic_analyzer.py`**: Retrieves the input FASTA filename from command-line arguments and calls the main function in `analyze_fasta.py`.
- **`analyze_fasta.py`**: Contains a main function and several helper functions that analyze FASTA sequences and print a summary for each sequence, including:
  - Sequence length
  - GC content
  - Reverse complement
  - Longest open reading frame (ORF) sequence[^1]

## 2. Revised Analyzer

The revised analyzer (`revised_analyzer.py` and `sample_qc.py`) builds upon the simple analyzer by introducing additional Python concepts and improving how FASTA records are processed.

Key improvements include:

- **Generator functions:** Streams FASTA records from the input file rather than loading the entire file into memory at once.
- **Object-oriented programming:** Organizes sequence analysis methods into a `Sample` class.
- **Comprehensions:** Uses list and other comprehensions to simplify data processing.
- **Error handling:** Wraps file reading and numeric parsing operations in `try/except` blocks so that malformed records can be handled without interrupting the processing of the remaining sequences.
- **Output file:** Saves the analysis summary to a results file rather than simply printing it to the terminal.

## 3. Testing

I used `pytest` in `test_analyzer.py` to verify that the revised analyzer correctly calculates GC content, generates reverse complements, and identifies the longest ORF for 10 sample sequences.

The test sequences were generated using `sample_fasta_generator.py`. Each sequence consists of 1,000 randomly generated nucleotide bases.

To manually establish the expected values for the tests, I used the following online tools:

| Tool | Purpose |
|---|---|
| [GC Content Calculator](https://jamiemcgowan.ie/bioinf/gc_content.html) | Calculates GC content |
| [Reverse Complement Tool](https://www.bioinformatics.org/sms/rev_comp.html) | Generates reverse complements |
| [NCBI ORF Finder](https://www.ncbi.nlm.nih.gov/orffinder/) | Identifies open reading frames |

[^1]: The longest ORF search algorithm in both the simple and revised analyzers only considers the forward strand. It does not search the reverse complement.
