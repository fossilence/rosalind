input_file = open('rosalind_orf.txt', 'r')
output_file = open('rosalind_orf_output', 'w')
fasta = input_file.read().replace('\n', '')
rna = fasta[14:].replace('T', 'U')
rev_rna = rna[::-1]

print(fasta + '\n')

print(rna + '\n')
print(rev_rna + '\n')

rc_rna = rev_rna.replace('A','u').replace('U','a').replace('G', 'c').replace('C', 'g').upper()

print(rc_rna +'\n')

trans_table = {'UUU' : 'F', 'CUU' : 'L', 'AUU' : 'I', 'GUU' : 'V', 'UUC' : 'F', 'CUC' : 'L', 'AUC' : 'I', 'GUC' : 'V', 'UUA' : 'L', 'CUA' : 'L', 'AUA' : 'I', 'GUA' : 'V', 'UUG' : 'L', 'CUG' : 'L', 'AUG' : 'M', 'GUG' : 'V', 'UCU' : 'S', 'CCU' : 'P', 'ACU' : 'T', 'GCU' : 'A', 'UCC' : 'S', 'CCC' : 'P', 'ACC' : 'T', 'GCC' : 'A', 'UCA' : 'S', 'CCA' : 'P', 'ACA' : 'T', 'GCA' : 'A', 'UCG' : 'S', 'CCG' : 'P', 'ACG' : 'T', 'GCG' : 'A', 'UAU' : 'Y', 'CAU' : 'H', 'AAU' : 'N', 'GAU' : 'D', 'UAC' : 'Y', 'CAC' : 'H', 'AAC' : 'N', 'GAC' : 'D', 'UAA' : 'x', 'CAA' : 'Q', 'AAA' : 'K', 'GAA' : 'E', 'UAG' : 'x', 'CAG' : 'Q', 'AAG' : 'K', 'GAG' : 'E', 'UGU' : 'C', 'CGU' : 'R', 'AGU' : 'S', 'GGU' : 'G', 'UGC' : 'C', 'CGC' : 'R', 'AGC' : 'S', 'GGC' : 'G', 'UGA' : 'x', 'CGA' : 'R', 'AGA' : 'R', 'GGA' : 'G', 'UGG' : 'W', 'CGG' : 'R', 'AGG' : 'R', 'GGG' : 'G'} 

rna_frames = []
rc_rna_frames = []

i = 0
while i < 3:
    rna_frames.append(rna[i:])
    rc_rna_frames.append(rc_rna[i:])
    i += 1

def split_codons(string):
    start, end = 0, 3
    codon_list = []
    while start < len(string) - 2:
        codon_list.append(string[start:end])
        start = start + 3
        end = end + 3
    return codon_list

all_frames = rna_frames + rc_rna_frames
all_codons = []

print(all_frames[5])

for seq in all_frames:
    all_codons.append(split_codons(seq))

for element in all_codons:
   for key, value in trans_table.items():
        i = 0
        while i < len(element):
            if key == element[i]:
             element[i] = value
            i += 1
print(all_codons[5])
aa_list = []
for i in all_codons:
    aa_list.append(''.join(i))

print('Problem frame: translated ')
print(aa_list[5])

protein_list = []
for element in aa_list:
    i = 0
    for char in element:
        if char == 'M':
            M = element.find('M', i)
            i += M + 1
            if element.find('x', M) != -1:
                x = element.find('x', M)
                protein_list.append(element[M:x])

print(protein_list)
answer = list(set(protein_list))
for i in answer:
    output_file.write(i + '\n')

input_file.close()
output_file.close()
