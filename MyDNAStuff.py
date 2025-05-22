#reversed codon and converts it to lowercase
def reverse_codon(codon):
    first = codon[0]
    second = codon[1]
    third = codon[2]
    reversed_codon = third+second+first
    return reversed_codon.upper()


#finds position and frame of first start codon (ATG)
def find_start_codon(dna_seq):
    atg_pos = dna_seq.lower().find('atg')
    frame = (atg_pos % 3) + 1 if atg_pos != -1 else None
    return (atg_pos+1,frame) if atg_pos != -1 else (None, None)


#analyzes gene sequences
def analyze_gene_seq(gene_seq):
    met_codon = 'ATG'
    starts_with_met = gene_seq.startswith(met_codon)
    met_pos = gene_seq.find(met_codon)
    met_frame = (met_pos % 3) + 1 if met_pos != -1 else None

    num_of_nucs = len(gene_seq)
    num_of_aa = (num_of_nucs // 3) - 1
    num_of_gcs = gene_seq.count('C') + gene_seq.count('G')
    gc_percent = 100 * (num_of_gcs/num_of_nucs)

    return {
        'Starts with MET:',starts_with_met,
        'MET Position:',met_pos + 1 if met_pos != -1 else None,
        'Met Frame:',met_frame,
        'Number of Nucleotides:',num_of_nucs,
        'Number of Amino Acids:',num_of_aa,
        'GC Percent:',round(gc_percent,2)
    }


#checks if sequence is tandem repeat
def tandem_repeat(seq):
    for chunk_len in range(1,len(seq)//2+1):
        chunk = seq[:chunk_len]
        copies = len(seq) // chunk_len
        if seq == chunk * copies:
            return (True, copies, chunk)
    return (False, 0, "")


#returns complement of nucleotide
def complement(nuc):
    complements = {'A':'T','T':'A','C':'G','G':'C','N':'N',
                   'a':'t','t':'a','c':'g','g':'c','n':'n'}
    return complements.get(nuc,nuc)


#returns reverse complement of DNA primer
def reverse_complement(seq):
    return ''.join(complement(n) for n in reversed(seq))


#checks if primer is reverse complement palindrome
def palindrome(primer):
    try:
        half_len = len(primer)//2
        first_half = primer[:half_len]
        second_half = primer[-half_len:]
        return first_half == reverse_complement(second_half)
    except TypeError:
        print('Input primer must be a string')
    except ValueError:
        print('Input primer cannot be empty')


#reads sequence from file
def readseq(seqfile):
    try:
        theseq = ''.join(open(seqfile).read().split())
        theseq = theseq.upper()
        return theseq
    except FileNotFoundError:
        print('Error: The file was not found')
    except IOError:
        print('Error reading file')


# read codons from a file
def readcodons(codonfile):
    f = open(codonfile)
    data = {}
    for l in f:
        sl = l.split()
        key = sl[0]
        value = sl[2]
        data[key] = value
    f.close()

    try:
        b1 = data['Base1']
        b2 = data['Base2']
        b3 = data['Base3']
        aa = data['AAs']
        st = data['Starts']
    except KeyError:
        print("Error: Missing required key in codon table")
        
    codons = {}
    init = {}
    n = len(aa)
    for i in range(n):
        codon = b1[i] + b2[i] + b3[i]
        codons[codon] = aa[i]
        init[codon] = (st[i] == 'M')
    return codons,init


#translates sequence using codon table
def trans(codons, seq):
    aalist = []
    for i in range(0,len(seq),3):
        try:
            codon = seq[i:i+3]
            aa = codons.get(codon,'X')
            aalist.append(aa)
        except IndexError:
            print('IndexError at position i.  Adding "X" for incomplete codon')
    aaseq = ''.join(aalist)
    return aaseq







