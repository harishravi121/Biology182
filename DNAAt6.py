import random
DNA='ATACAG'
w='ATCG'
s=''
N=random.randint(5,10)
for j in range(N):
    a=random.sample(w,1)
    s=s+a[0]
DNA=s
print(DNA)
A=[5,5,5,0] #CHNO
T=[5,6,2,2]
C=[4,5,3,1]
G=[5,5,5,1]
DNAf=[15,31,3,13,2]
x1=[0,0,0,0]
def formula1(DNA1):
    global x1;
    An=0
    Tn=0
    Cn=0
    Gn=0
    for i in range(len(DNA1)):
        if(DNA1[i]=='A'):
            An=An+1
        if(DNA1[i]=='T'):
            Tn=Tn+1
        if(DNA1[i]=='C'):
            Cn=Cn+1
        if(DNA1[i]=='G'):
            Gn=Gn+1

    x1=[An, Tn ,Cn,Gn,len(DNA1)]


#print(x1)
def Molecules():
    global A,T,C,G,DNAf
    #print(A)
    C1=x1[0]*A[0]+x1[1]*T[0]+x1[2]*C[0]+x1[3]*G[0]+x1[4]*DNAf[0]
    H=x1[0]*A[1]+x1[1]*T[1]+x1[2]*C[1]+x1[3]*G[1]+x1[4]*DNAf[1]
    N=x1[0]*A[2]+x1[1]*T[2]+x1[2]*C[2]+x1[3]*G[2]+x1[4]*DNAf[2]
    O=x1[0]*A[3]+x1[1]*T[3]+x1[2]*C[3]+x1[3]*G[3]+x1[4]*DNAf[3]
    P=x1[4]*DNAf[4]
    HC='C'+str(C1)+'H'+str(H)+'N'+str(N)+'-O_'+str(O)+'P'+str(P)
    print('DNA Molecule  ', HC)
flag=0
NAME="HIR VPROF DR HARISH RAVI"
#print('Bio Lesson')
#print('Proteins make amino acids.There are about 20 Amino Acids, three base pairs code for an amino acid which can be represented with A-Z and 6 extra dummy acids. An extra base pair X has to be added to code for some letters as some codes are redundant and blanks are hard to find. Eg : the triplet ATA codes for I, there is ATCG base pairs and complementary in where A=T and C=G double strand. RNA has U instead of T.')
#print('There is 3D folding and active sites which is great')

w='ATCG'
flag=0
def translate(seq): 
       
    table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
    } 
    protein ="" 
    if len(seq)%3 == 0: 
        for i in range(0, len(seq), 3): 
            codon = seq[i:i + 3] 
            protein+= table[codon] 
    return protein 

def revtranslate(seq): 
       
    table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
        'XXX':' ', 'ATX':'U','AXT':'O','XAG':'X','XTG':'B',
        'XTA':'J','XTA':'J','XGA':'Z','.':'.'
    } 
    DNA ="" 
    table2 = {v: k for k, v in table.items()}
    for i in range(0, len(seq)): 
         
        DNA+= table2[seq[i]] 
    return DNA
def RNA(seq): 
       
    table3 = { 
        'A':'U','T':'A','C':'G','G':'C','X':'X','.':'.'
    } 
    RNA ="" 
    
    for i in range(0, len(seq)): 
         
        RNA+= table3[seq[i]] 
    return RNA


protein=NAME
print('Protein ',protein,'length ', len(protein))
DNA=revtranslate(protein)
RNA=RNA(DNA)
print('DNA ',DNA, 'length ', len(DNA))

print('RNA ',RNA,'length ', len(RNA))

formula1(DNA)
Molecules()


