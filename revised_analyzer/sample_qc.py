from typing import Final

STOP_CODONS: Final = ("TAA", "TAG", "TGA")

class Sample:
    def __init__(self, header, sequence):
        self.sample_id = header
        self.sample_seq = sequence

    def get_gc_content(self):
        #try/except to handle divided by zero error if the sequence is empty
       try:
           gc_count = self.sample_seq.count("G") + self.sample_seq.count("C")
           total = len(self.sample_seq)
           gc_content = (gc_count / total) * 100.0
       except ZeroDivisionError:
           print(f"Error: Sequence for sample '{self.sample_id}' is empty. Cannot calculate GC content.")
           return 0.0
       return gc_content 

    def get_reverse_complement(self):
        #use try/except to handle invalid input for the sequence
        try:
            complement_map = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
            complemented = "".join(complement_map[base] for base in self.sample_seq)
            return complemented[::-1]
        except KeyError as e:
            print(f"Error: Invalid base '{e.args[0]}' found in sequence for sample '{self.sample_id}'. Cannot compute reverse complement.") 
            return "Error: Invalid sequence"

    def get_longest_ORF_seq(self):
        #TODO: Handle equal length?
        longest_ORF_seq = ""
        for i in range(3):
            curr_orf_seq = self.sample_seq[i:]
            curr_orf_longest_seq = self.__get_longest_seq(curr_orf_seq)
            if len(curr_orf_longest_seq) > len(longest_ORF_seq):
                longest_ORF_seq = curr_orf_longest_seq
        if len(longest_ORF_seq) == 0:
            longest_ORF_seq = "No sequences found in any ORF"
        return longest_ORF_seq

    def __get_longest_seq(self, sequence):
        start_indices = [i for i in range(len(sequence)) if sequence.startswith("ATG", i)]
        if not start_indices:
            return ""
        longest_seq = ""
        for start_index in start_indices:
            curr_seq = sequence[start_index:]
            for i in range(start_index, len(curr_seq) - 2, 3):
                curr_codon = curr_seq[i:i+3]
                if (curr_codon in STOP_CODONS) and (len(curr_seq[0: i + 3]) > len(longest_seq)):
                    longest_seq = curr_seq[0:i + 3]
                    break
        return longest_seq
        