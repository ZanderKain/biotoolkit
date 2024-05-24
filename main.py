from DNATOOLKit import *
from utilities import colored 
from structures import *

import random

# Create a random DNA sequence 
rndDNAstr = ''.join([random.choice(nucleotides) for nuc in range(50)])

dnaStr = validateSeq(rndDNAstr)

print(f'Sequence: {(dnaStr)}\n')
print(f'[1] + Sequence Length: {len(dnaStr)}\n')
print((f'[2] + Nucleotide Frequency: {countNucFrequency(dnaStr)}\n'))
print(f'[3] + Transcription of DNA to RNA: {(transcription(dnaStr))}\n')

print(f"[4] + DNA String + Reverse Complement:\n5' {(dnaStr)} 3' ")
print(f"   {''.join(['|'for counts in range(len(dnaStr))])}")
print(f"5' {reverse_complement(dnaStr)[::-1]} 3' *Complement*")
print(f"3' {reverse_complement(dnaStr)} 5' *Reverse Complement*")

print(f"[5] + GC contents in percentage of sequenced string: {gc_content(dnaStr)}%")

print(f"[6] + Subsequence of GC content in k = 10: {gc_subseq_content(dnaStr, k=10)}")

print(f"[7] + Aminoacids Sequence from DNA: {translate_seq(dnaStr)}\n")

print(f"[8] + Codon frequency (L): {codon_usage(dnaStr, 'L')}\n")

print(f'[9] + reading Frames:')
for frame in gen_reading_frames(dnaStr):
    print(frame)

print(f'\n[10] + All Proteins in 6 open reading frames:') 
for prot in all_proteins_from_orfs(dnaStr,0,0,True):
    print(f'{prot}')