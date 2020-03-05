file = open('rosalind_orf.txt', 'r')
fasta = file.readlines()
rna = fasta[1].replace('T', 'U')
rev_rna = fasta[1][::-1].replace('T', 'U')


rc_rna = rev_rna.replace('A','u').replace('U','a').replace('G', 'c').replace('C', 'g').upper()


print('DNA')
print(fasta[1])

print('RNA')
print(rna)
print(rna[0:])
print(rna[1:])
print(rna[2:])
print('Reverse complement RNA')
print(rc_rna)
print(rc_rna[0:])
print(rc_rna[1:])
print(rc_rna[2:])

trans_table = {'UUU' : 'F', 'CUU' : 'L', 'AUU' : 'I', 'GUU' : 'V', 'UUC' : 'F', 'CUC' : 'L', 'AUC' : 'I', 'GUC' : 'V', 'UUA' : 'L', 'CUA' : 'L', 'AUA' : 'I', 'GUA' : 'V', 'UUG' : 'L', 'CUG' : 'L', 'AUG' : 'M', 'GUG' : 'V', 'UCU' : 'S', 'CCU' : 'P', 'ACU' : 'T', 'GCU' : 'A', 'UCC' : 'S', 'CCC' : 'P', 'ACC' : 'T', 'GCC' : 'A', 'UCA' : 'S', 'CCA' : 'P', 'ACA' : 'T', 'GCA' : 'A', 'UCG' : 'S', 'CCG' : 'P', 'ACG' : 'T', 'GCG' : 'A', 'UAU' : 'Y', 'CAU' : 'H', 'AAU' : 'N', 'GAU' : 'D', 'UAC' : 'Y', 'CAC' : 'H', 'AAC' : 'N', 'GAC' : 'D', 'UAA' : 'x', 'CAA' : 'Q', 'AAA' : 'K', 'GAA' : 'E', 'UAG' : 'x', 'CAG' : 'Q', 'AAG' : 'K', 'GAG' : 'E', 'UGU' : 'C', 'CGU' : 'R', 'AGU' : 'S', 'GGU' : 'G', 'UGC' : 'C', 'CGC' : 'R', 'AGC' : 'S', 'GGC' : 'G', 'UGA' : 'x', 'CGA' : 'R', 'AGA' : 'R', 'GGA' : 'G', 'UGG' : 'W', 'CGG' : 'R', 'AGG' : 'R', 'GGG' : 'G'} 

rna_frames = []
rc_rna_frames = []

print(len(trans_table))

i = 0
while i < 3:
    rna_frames.append(rna[i:])
    i += 1

i = 0
while i < 3:
    rc_rna_frames.append(rc_rna[i:])
    i += 1

print('RNA frames')
for element in rna_frames:
    print(element)
#    print('Codons: ' + str(len(element)/3))

print('RC RNA frames')
for element in rc_rna_frames:
    print(element)
#    print('Codons: ' + str(len(element/3))

def split_codons(string):
    start, end = 0, 3
    codon_list = []
    while start < len(string) - 2:
        codon_list.append(string[start:end])
        start = start + 3
        end = end + 3
    return codon_list

c = [] #all 6 reading frames in a list of lists
i = 0
while i < 3:
    c.append(split_codons(rna_frames[i]))
    i += 1

i = 0
while i< 3:
    c.append(split_codons(rc_rna_frames[i]))
    i += 1

print(c)

print('List of all reading frames')
for i in c:
    print(len(i))

for element in c:
   for key, value in trans_table.items():
        i = 0
        while i < len(element):
            if key == element[i]:
             element[i] = value
            i += 1
            
# c contains a list of all translated amino acids into a list
#time to search for a 'M' and end on a 'Stop' and return the all the aminos inbetween
#turn list of lists into list of strings

aa_list = []
i = 0
while i < len(c):
    aa_list.append(''.join(c[i]))
    i += 1

print(aa_list)
protein_list = []
for element in aa_list:
    i = 0
    for char in element:
        if char == 'M':
            M = element.find('M', i)
            i = i + M  + 1
            if element.find('x', i) != -1:
                x = element.find('x')
                protein_list.append(element[M:x])

print(protein_list)
#for element in protein_list:
#    if element[0] == 'M':
#        print(element)


