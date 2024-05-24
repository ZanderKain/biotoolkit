# DNA Toolkit File
from structures import *
from collections import Counter
# Check the sequence to make sure it is a DNA String
def validateSeq(dna_seq):
    tempseq = dna_seq.upper()
    for nuc in tempseq:
        if nuc not in nucleotides:
            return False
    return tempseq

# Counts the Nucleotides in a DNA String
def countNucFrequency(seq):
    temp_freq_dict = {"A": 0, "C": 0, "G": 0, "T": 0}
    for nuc in seq:
        temp_freq_dict[nuc] += 1
    return temp_freq_dict

# Transcribes DNA into RNA (T->U)
def transcription(seq):
    return seq.replace("T","U")

def reverse_complement(seq):
    # return "".join([DNA_reverse_complement[nuc] for nuc in seq])
    #Second way to solve the reverse complement. Faster
    
    translated_seq = str.maketrans("ATCG", "TAGC")
    return seq.translate(translated_seq)[::-1]

def gc_content(seq):
    return round((seq.count('C') + seq.count('G')) / len(seq) * 100)

def gc_subseq_content(seq, k=20):
    # GC content in a DNA/RNA sub-sequence length k, k=20 by default.
    gc_list = []
    # Cuts the DNA sequence from 0 index to last index at a step of k
    for i in range(0,len(seq) - k + 1, k):
        # Loops through sub sequence
        subseq = seq[i:i + k]
        gc_list.append(gc_content(subseq))
    return gc_list

# Codon Translater 
def translate_seq(seq, init_pos=0):
    """Translates a DNA sequence into an aminoacid sequence"""
    return [DNA_codons[seq[pos:pos +3]] for pos in range(init_pos,len(seq)-2, 3)]

def codon_usage(seq, aminoacid):
    """Provides the frequency of each codon encoding a given aminoacid in a DNA sequence"""
    tmplist = []
    for i in range(0, len(seq)-2,3):
        if DNA_codons[seq[i:i+3]] == aminoacid:
            tmplist.append(seq[i:i+3])
    
    freqdict = dict(Counter(tmplist))
    totalwight = sum(freqdict.values())
    for seq in freqdict:
        freqdict[seq] = round(freqdict[seq]/ totalwight, 2)
    return freqdict

def gen_reading_frames(seq):
    """Generate the six reading frames of a DNA sequence, including the reverse complement"""
    frames = []
    frames.append(translate_seq(seq,0))    
    frames.append(translate_seq(seq,1))    
    frames.append(translate_seq(seq,2))    
    frames.append(translate_seq(reverse_complement(seq),0))    
    frames.append(translate_seq(reverse_complement(seq),1))    
    frames.append(translate_seq(reverse_complement(seq),2))    

    return frames

def proteins_from_rf(aa_seq):
    """Compute all possible proteins in an aminoacid seq and return a list of possible proteins"""
    current_prot =[]
    proteins = []
    
    for aa in aa_seq:
        if aa == "_":
            #Stop accumlating amino acid
            if current_prot:
                for p in current_prot:
                    proteins.append(p)
                current_prot =[] 
        
        else:
            # Start accumlating amino acids
            if aa == "M":
                current_prot.append("")
            for i in range(len(current_prot)):
                current_prot[i] += aa
    return proteins
   
def all_proteins_from_orfs(seq, startReadPos=0, endReadPos=0,ordered=False):
    """Compute all possible proteins for all open reading frames"""
    """Protein Search Database: https://www.ncbi.nlm.nih.gov/nuccore/nM_00118509"""
    """API can be used to pull protein info""" 
    if endReadPos > startReadPos:
        rfs = gen_reading_frames(seq[startReadPos:endReadPos])
    else:
        rfs = gen_reading_frames(seq)
    
    res = []
    for rf in rfs:
        prots = proteins_from_rf(rf)
        for p in prots:
            res.append(p)
    
    if ordered:
        return sorted(res, key=len, reverse=True)