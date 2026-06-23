from Bio import SeqIO
import sys

with open(sys.argv[2], "w") as out:
    for record in SeqIO.parse(sys.argv[1], "fasta"):
        n_count = str(record.seq).upper().count("N")
        if n_count <= 50:
            SeqIO.write(record, out, "fasta")
