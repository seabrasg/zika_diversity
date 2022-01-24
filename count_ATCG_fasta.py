# This script calculates, for each sequence, the percentage of sites ATCG (not ambiguous) in the total coding region (10269 bp)
# It also calculates GC content

# usage: python3 count_ATCG_fasta.py input.fasta output.txt

import sys
from Bio import SeqIO

fasta_file = sys.argv[1]  # Input fasta file with all records
result_file = sys.argv[2] # Output txt file

wanted = set()

fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')

with open(result_file, "w") as f:
    f.write('Accession\tATCG%\tCG%\n')

    for r in fasta_sequences:
        accession=r.id.split(".")[0]
        A_count=r.seq.count('A')
        C_count=r.seq.count('C')
        G_count=r.seq.count('G')
        T_count=r.seq.count('T')
        length=len(r.seq)
        ATCG_percentage=float(A_count+C_count+G_count+T_count)*100.0/10269
        CG_percentage=float(C_count+G_count)*100.0/(A_count+C_count+G_count+T_count)
        output_line = '%s\t%f\t%f\n' % (accession, ATCG_percentage, CG_percentage)

        f.write(output_line)
    f.close()
 
