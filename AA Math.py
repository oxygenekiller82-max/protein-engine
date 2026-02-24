#Formulas file 
import math 

def calculate_hydrophobicity(sequence):
    """ 1/n * Sum(Hi): returns hydrophobicity result of the list of AA"""
    if len(sequence)==0: return 0 
    total_hydro=sum(aa.hydro for aa in sequence)
    return total_hydro/len(sequence)


def calcualte_entropy(seqeuence):
    """-Sum (pi*log2(pi)), i = UNIQUE amino acid in the SET"""
    if len(seqeuence)==0: return 0

    entropy=0

    set_count={} #UNIQUE AAs + their count 
    for aa in seqeuence:
        set_count[aa.code] = set_count.get(aa.code,0)+1

    for code in set_count:
        p_i= set_count[code]/len(seqeuence) #FREQUENCE de chaque AA UNQIUE
        entropy-=p_i*math.log2(p_i)

    return entropy



