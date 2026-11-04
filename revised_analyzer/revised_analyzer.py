import sample_qc

s1 = sample_qc.Sample("header", "ATGCCATGGCATTAGCCATGGCATAAACATGGCTGA")
orf = s1.get_longest_ORF_seq()
print(orf)