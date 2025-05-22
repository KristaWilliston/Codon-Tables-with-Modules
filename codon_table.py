import sys

#reads DNA seq from fiiel and returns it as uppercase string
def read_seq_from_filename(seqfile):
    try:
        theseq = ''.join(open(seqfile).read().split())
        theseq = theseq.upper()
        return theseq
    except FileNotFoundError:
        print("Error: The file was not found.")
        sys.exit(1)
    except IOError:
        print("Error reading file")
        sys.exit(1)


#reads codon table from file and returns dictionary of codons and their amino acids
def read_codons_from_filename(codonfile):
    f = open(codonfile)
    data = {}
    for line in f:
        sl = line.split()
        key = sl[0]
        value = sl[2]
        data[key] = value
    f.close()

    b1 = data['Base1']
    b2 = data['Base2']
    b3 = data['Base3']
    aa = data['AAs']
    st = data['Starts']

    codons = {}
    init = {}
    for i in range(len(aa)):
        codon = b1[i] + b2[i] + b3[i]
        codons[codon] = aa[i]
        init[codon] = (st[i] == 'M')
    return codons, init


#checks if given codon is an initiation codon
def is_init(codons, seq):
    try:
        return seq[:3] in codons and codons[seq[:3]]
    except TypeError:
        print('Input sequence must be a string')
    except ValueError:
        print('Input sequence must be at least 3 characters long')

 
#reverse complement of DNA seq
def reverse_complement(seq):
    comp_dict = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
    try:
        return ''.join(comp_dict.get(nuc,nuc) for nuc in reversed(seq))
    except KeyError:
        print('Error: Invalid nucleotide in sequence')


#returns amino acid fro codon with ambiguity or 'X' if multiple aa's possible
def get_ambig_aa(codon_table, codon):
    aas = set()
    for n3 in 'ACGT':
        codon1 = codon[:2] + n3
        aas.add(codon_table.get(codon1, 'X'))
    return 'X' if len(aas) > 1 else aas.pop().upper()


#translate seq using codon table, accommodating for frames and ambiguity
def translate(codons, seq, frame):
    seq += 'N' * (frame-1)
    aa_list = []

    for i in range(frame-1, len(seq), 3):
        codon = seq[i:i+3]

        if codon in codons:
            aa = codons[codon]
        elif codon.count('N') == 1 and codon[2] == 'N':
            aa = get_ambig_aa(codons,codon)
        else:
             aa = 'X'

        aa_list.append(aa)
    aaseq = ''.join(aa_list)
    return aaseq


def translate_all_frames(codons, seq):
    translations = {}
    
    for frame in range(1,4):
        translations[frame] = translate(codons,seq,frame)

    rev_seq = reverse_complement(seq)

    for frame in range(1,4):
        translations[frame + 3] = translate(codons,rev_seq,frame)

    return translations
