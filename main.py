from Bio import SeqIO

filekk = SeqIO.parse("dna.example.fasta", "fasta")

records = [rec for rec in filekk]

print(len(records))



# for seq_record in SeqIO.parse("dna.example.fasta", "fasta"):
#     print(seq_record.id)
#     print(repr(seq_record.seq))
#     print(len(seq_record))