import sys
import os
import re

#filename = sys.argv[1]
seq = 'atgaaaaagaatacattaagtgcaatattaatgactttatttttatttatatcttgtaataattcagggaaagatgggaatacatctgcaaattctgctgatgagtctgttaaagggcctaatcttacagaaataagtaaaaaaattacggattctaatgcggttttacttgctgtgaaagaggttgaagcgttgctgtcatctatagatgaaattgctgctaaagctattggtaaaaaaatacaccaaaataatggtttggataccgaaaataatcacaatggatcattgttagcgggagcttatgcaatatcaaccctaataaaacaaaaattagatggattgaaaaatgaaggattaaaggaaaaaattgatgcggctaagaaatgttctgaaacatttactaataaattaaaagaaaaacacacagatcttggtaaagaaggtgttactgatgctgatgcaaaagaagccattttaaaaacaaatggtactaaaactaaaggtgctgaagaacttggaaaattatttgaatcagtagaggtcttgtcaaaagcagctaaagagatgcttgctaattcagttaaagagcttacaagccctgttgtggcagaaagtccaaaaaaaccttaa'

AAtable={
        'TTT':'Phe','TTC':'Phe','TTA':'Leu','TTG':'Leu',
        'CTT':'Leu','CTC':'Leu','CTA':'Leu','CTG':'Leu',
        'ATT':'Ile','ATC':'Ile','ATA':'Ile','ATG':'Met',
        'GTT':'Val','GTC':'Val','GTA':'Val','GTG':'Val',
        'TCT':'Ser','TCC':'Ser','TCA':'Ser','TCG':'Ser',
        'CCT':'Pro','CCC':'Pro','CCA':'Pro','CCG':'Pro',
        'ACT':'Thr','ACC':'Thr','ACA':'Thr','ACG':'Thr',
        'GCT':'Ala','GCC':'Ala','GCA':'Ala','GCG':'Ala',
        'TAT':'Tyr','TAC':'Tyr','TAA':'***','TAG':'***',
        'CAT':'His','CAC':'His','CAA':'Gln','CAG':'Gln',
        'AAT':'Asn','AAC':'Asn','AAA':'Lys','AAG':'Lys',
        'GAT':'Asp','GAC':'Asp','GAA':'Glu','GAG':'Glu',
        'TGT':'Cys','TGC':'Cys','TGA':'***','TGG':'Trp',
        'CGT':'Arg','CGC':'Arg','CGA':'Arg','CGG':'Arg',
        'AGT':'Ser','AGC':'Ser','AGA':'Arg','AGG':'Arg',
        'GGT':'Gly','GGC':'Gly','GGA':'Gly','GGG':'Gly'}
#print(AAtable)

non_count = 0
syn_count = 0

seq = seq.upper()
codonnum = len(seq)//3
index = codonnum*3

for i in range(0, index, 3):
        codon = seq[i:i + 3]
        protein = AAtable[codon]
        #print(protein)
        for i in range(3):
                if codon[i] != 'A':
                        new_codon = codon.replace(codon[i], 'A')
                        new_protein = AAtable[new_codon]
                        if new_protein == protein:
                                syn_count += 1
                        else:
                                non_count += 1
                if codon[i] != 'C':
                        new_codon = codon.replace(codon[i], 'C')
                        new_protein = AAtable[new_codon]
                        if new_protein == protein:
                                syn_count += 1
                        else:
                                non_count += 1
                if codon[i] != 'T':
                        new_codon = codon.replace(codon[i], 'T')
                        new_protein = AAtable[new_codon]
                        if new_protein == protein:
                                syn_count += 1
                        else:
                                non_count += 1
                if codon[i] != 'G':
                        new_codon = codon.replace(codon[i], 'G')
                        new_protein = AAtable[new_codon]
                        if new_protein == protein:
                                syn_count += 1
                        else:
                                non_count += 1
print(syn_count)
print(non_count)








