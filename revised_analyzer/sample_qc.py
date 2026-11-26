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

    # Comment out the existing print statements for debugging purposes
    def get_longest_ORF_seq(self):
        #TODO: Handle equal length?
        longest_ORF_seq = ""
        for frame in range(3):
            # print(f"CHECKING FRAME {frame}---------------------------------------------")
            curr_orf_longest_seq = self.__get_longest_seq(self.sample_seq, frame)
            if len(curr_orf_longest_seq) > len(longest_ORF_seq):
                longest_ORF_seq = curr_orf_longest_seq
                # print(f"New longest ORF found in frame {frame}")
                # print(f"CURRENT ORF seq: {longest_ORF_seq}")
        if len(longest_ORF_seq) == 0:
            longest_ORF_seq = "No sequences found in any ORF"
        return longest_ORF_seq
# Comment out the existing print statements for debugging purposes
    def __get_longest_seq(self, sequence, frame):
        start_indices = [i for i in range(frame, len(sequence) - 2, 3) if sequence.startswith("ATG", i)]
        if not start_indices:
            return "No start codons found in this ORF"
        longest_seq = ""
        best_start_index = -1
        for start_index in start_indices:
            # print(f"Start index: {start_index}")
            for i in range(start_index, len(sequence) - 2, 3):
                curr_codon = sequence[i:i+3]
                if (curr_codon in STOP_CODONS):
                    # print(f"STOP codon: {curr_codon}")
                    curr_seq = sequence[start_index:i+3]
                    if len(curr_seq) > len(longest_seq):
                        longest_seq = curr_seq
                        best_start_index = start_index
                    break
        # print(f"Best start index: {best_start_index}")
        # print(f"Longest OF THIS FRAME: {longest_seq}")
        return longest_seq

        