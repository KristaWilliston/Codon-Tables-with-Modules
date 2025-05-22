# Codon-Tables-with-Modules

The result is the translation of all 6 frames of a DNA sequence by utilizing my MyDNAStuff and codon_table modules, my main script ModuleForCodonTables, and a DNA sequence from the commandline.

For both MyDNAStuff and codon_tables modules, the solutions added were modified so that everything in them all became defined functions and were ambiguous enough that they would work with any string of DNA sequence.  I had to add to my codon_table module a new translate function that would translate the input DNA sequence in 6 frames.  I then created a main script that modified the code to read from each of these modules by calling them from the command line and printed out the frame and amino acid sequence of the input DNA sequence.

