from MyDNAStuff import *
from codon_table import *
import sys

if len(sys.argv) < 3:
    print("Require codon table and DNA sequence on command-line.")
    sys.exit(1)

table, _ = read_codons_from_filename(sys.argv[1])
seq = read_seq_from_filename(sys.argv[2])

if is_init(table,seq):
    print("Initial codon is a translation start codon")
#print(table)
#print(type(table))
translations = translate_all_frames(table,seq)

for frame, seq in translations.items():
    print("Frame",frame,":",seq)
