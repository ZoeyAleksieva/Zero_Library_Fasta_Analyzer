from collections import defaultdict

def analyze_fasta():
    #key : header, value = [length, GC content,reverse complement, longest open reading frame (ORF)]
    #default dicts automatically create empty list when a new key is accessed

    sequence_info = defaultdict(list)

    #Goal: read through fasta file and add to dictionary {header : [sequence]}
    #Start going through file and look for >
    #If >, then that's the start of a new sequence

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
    #Read by 3s and put in list
    codons = []

    #ONLY FOR 0TH FRAME, refactor to do 1st and 2nd
    #the -2 makes it stop before it reaches char that would make codon less than 3 chars
    for i in range(0, len(sequence) - 2, 3): #if len is 12, i is 0,3,6,9
        codons.append(sequence[i:i+3])
    
    atg_indices = []
    for i, codon in enumerate(codons):
        if codon == "ATG":
            atg_indices.append(i)
    print(f"ATG indices: {atg_indices}")
    stop_codons = {"TAA", "TAG", "TGA"}
    max_index_diff = 0
    curr_longest = []
    for index in atg_indices:
        curr_codons = codons[index:]
        #Once the first stop is found, you can break
        for i, codon in enumerate(curr_codons):
            if codon in stop_codons:
                print(f"found match: {codon}")
                if i - index > max_index_diff:
                    print(f"found new max: {index} to {i + index}")
                    max_index_diff = i - index
                    print(f"curr max: {codons[index : i + index + 1]}")
                    curr_longest = codons[index : i + index + 1]
                break
    #longest_ORF_codons = codons[max_atg_index : max_atg_index + max_index_diff]
    longest_ORF_seq = "".join(curr_longest)
    return longest_ORF_seq

analyze_fasta()