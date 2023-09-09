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

formula1(DNA)
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
Molecules()


