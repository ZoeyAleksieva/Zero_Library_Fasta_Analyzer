from collections import defaultdict
from typing import Final

STOP_CODONS: Final = ("TAA", "TAG", "TGA")

def analyze_fasta():
    
    sequence_info = defaultdict(list)

    with open("orf_test.fa", "r") as fasta_file:
        
        for line in fasta_file:
            line = line.strip()

            if line.startswith(">"):
                header = line[1:]
            else:
                if len(sequence_info[header]) == 0:
                    sequence_info[header].append(line)
                else: 
                    sequence_info[header][0] += line

        analyze_sequences(sequence_info)
        print_summary(sequence_info)
    return 

def analyze_sequences(sequence_info):
    for seq_id, info in sequence_info.items():
        sequence_info[seq_id].append(len(info[0]))
        sequence_info[seq_id].append(get_gc_content(info[0]))
        sequence_info[seq_id].append(get_reverse_complement(info[0]))
        sequence_info[seq_id].append(get_longest_ORF(info[0]))

def get_gc_content(sequence):
        gc_count = sequence.count("G") + sequence.count("C")
        gc_content = f"{(gc_count / len(sequence)) * 100:.2f}"
        return gc_content

def get_reverse_complement(sequence):
    complement_map = {"A": "T", "T": "A", "G": "C", "C": "G"}
    complemented = "".join(complement_map[base] for base in sequence)
    return complemented[::-1]

def print_summary(sequence_info):
    for seq, info in sequence_info.items():
        print("ID: ")
        print(seq)
        print("Info: ")
        print(info)

def get_longest_ORF(sequence):
    #call function for each ORF and then compare returned
    longest_ORF_seq = ""
    for i in range(3):
        print(f"Current ORF: {i}")
        codons = get_codons_for_ORF(i, sequence)
        curr_ORF_longest_seq = get_longest_seq(codons)
        print(f"Longest seq of curr ORF: {curr_ORF_longest_seq}")
        if len(curr_ORF_longest_seq) > len(longest_ORF_seq):
            longest_ORF_seq = curr_ORF_longest_seq
    return longest_ORF_seq

def get_codons_for_ORF(start_index, sequence):
    codons = []
    for i in range(start_index, len(sequence) - 2, 3): #if len is 12, i is 0,3,6,9
        codons.append(sequence[i:i+3])
    return codons

def get_longest_seq(codons):
    atg_indices = get_ATG_indices(codons)
    max_index_diff = 0
    curr_longest = []
    for index in atg_indices:
        curr_codons = codons[index:]
        for i, codon in enumerate(curr_codons):
            if codon in STOP_CODONS:
                print(f"found match: {codon}")
                if i - index > max_index_diff:
                    print(f"found new max: {index} to {i + index}")
                    max_index_diff = i - index
                    print(f"curr max: {codons[index : i + index + 1]}")
                    curr_longest = codons[index : i + index + 1]
                break
    longest_seq = "".join(curr_longest)
    return longest_seq

def get_ATG_indices(codons):
    atg_indices = []
    for i, codon in enumerate(codons):
        if codon == "ATG":
            atg_indices.append(i)
    return atg_indices

analyze_fasta()