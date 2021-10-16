from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import matplotlib.pyplot as plt
def getmutationhotspot(protein_length,file,protein_num):
    b=[]
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
        b.append(result[0]/protein_num)
    plt.plot(range(protein_length),b)
        