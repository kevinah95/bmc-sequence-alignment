import numpy
seq1=["G","A","T","C","C","A"]
seq2=["G","T","G","C","C","T"]
gap=-2
match=1
mismatch=-1

matriz = numpy.zeros((5, 5))

first=0
for i in range(5):
    matriz[i][0]=first
    first=first+gap
first=0
for j in range(5):
    matriz[0][j]=first
    first=first+gap
    

def check_diago(value,posi,posj,seq1,seq2,match,mismatch):
    if(seq2[posi-1]==seq1[posj-1]):
        value=value+match
    else:
        value=value+mismatch
    return value


def check_up_left(value,gap):
    value=value+gap
    return value


for i in range(1,5):
    for j in range(1,5):
        diago=check_diago(matriz[i-1][j-1],i,j,seq1,seq2,match,mismatch)
        up=check_up_left(matriz[i-1][j],gap)
        left=check_up_left(matriz[i][j-1],gap)
        matriz[i][j]=max(diago,up,left)
        
print(matriz)