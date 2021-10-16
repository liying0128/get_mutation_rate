from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import matplotlib.pyplot as plt
def getmutationhotspot(file):
    similarity=[]
    Protein_length=[]
    paralog_num=0
    for i in SeqIO.parse(file,'fasta'):
        paralog_num=paralog_num+1
        Protein_length.append(len(i))
    protein_length=max(Protein_length)
    for i in range(protein_length):
        a=[]
        for seq_record in SeqIO.parse(file,'fasta'):
            seq=str(seq_record.seq)
            seq=seq.replace('*','')
            seq=Seq(seq)
            if len(seq)<i+1:
                a.append(0)
            else:
                a.append(seq[i])
        result=pd.value_counts(a)
        similarity.append((1-result[0]/paralog_num)*100)
    plt.plot(range(protein_length),similarity,color='red')
    plt.ylabel('Mutation rate (%)',fontsize=12)
    plt.xlabel('Position',fontsize=12) 