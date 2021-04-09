from Bio import SeqIO
from Bio.Seq import Seq
fasta_file = "/home/muggli/git2/Assemblies-of-putative-SARS-CoV2-spike-encoding-mRNA-sequences-for-vaccines-BNT-162b2-and-mRNA-1273/Figure1Figure2_032321.fasta"

# stolen from https://stackoverflow.com/questions/48799955/find-the-hamming-distance-between-two-dna-strings
def hamming_distance(s1, s2):
        if len(s1) != len(s2):
                    raise ValueError("Strand lengths are not equal!")
        return sum(ch1 != ch2 for ch1,ch2 in zip(s1,s2)) # if ch1 != "*" and ch2 != "*")
        # else:
        #     print("different lengths error in hamming()", len(s1),len(s2))

for seqno, seq_record in enumerate(SeqIO.parse(fasta_file, "fasta")):
#    print(dir(seq_record.seq))
    offset = seq_record.seq.find(Seq("ATG")) # start codon
    offset2 = seq_record.seq.find(Seq("TGATGA")) # stop codon
    offset2 += 6
#    offset2 = seq_record.seq.find(Seq("TGATGA")) # stop codon
    print("Seq:", seq_record.id)
    print("Start codon at:", offset)
    print("Stop codon at:", offset2)
    if seqno == 2:
        offset=0
        offset2 = len(seq_record.seq)
    seq = seq_record.seq[offset:offset2]
    print("Translating:", seq[:18],"...",seq[-18:])
    tr = seq.translate(table=2)
    print("done")
    print(tr)
    if seqno >= 1:
        print("translated diffs from first:", hamming_distance(seq.translate(table=2), last_seq.translate(table=2)), "out of", len(seq.translate(table=2)), "amino acids.")
        print("untranslated diffs from first:", hamming_distance(seq, last_seq), "out of",
              len(seq), "bases.")
    last_seq = seq

