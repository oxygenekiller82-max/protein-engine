#Formulas file 
import math 
from AA_BD import LIMITS
from AA_BD import BULKY_AAS

def calculate_hydrophobicity(sequence):
    """ 1/n * Sum(Hi): returns hydrophobicity result of the list of AA"""
    if len(sequence)==0: return 0 
    total_hydro=sum(aa.hydro for aa in sequence)
    return total_hydro/len(sequence)


def calcualte_entropy(sequence):
    """-Sum (pi*log2(pi)), i = UNIQUE amino acid in the SET"""
    if len(sequence)==0: return 0

    entropy=0

    set_count={} #UNIQUE AAs + their count 
    for aa in sequence:
        set_count[aa.code] = set_count.get(aa.code,0)+1

    for code in set_count:
        p_i= set_count[code]/len(sequence) #FREQUENCE de chaque AA UNQIUE
        entropy-=p_i*math.log2(p_i)

    return entropy

#Rq: Generallement un Indice Aliphatique >70 est considéré stable, <40 est non resistant 
#à l'augmentation de la temperature

#50-70 = "Mesophilic": Stable at body temperature (37C)
#>85 = Thermophilic: can survive higher tempratrures 
#<40 = Unstable: folds at even room temperture
def calculate_stability(sequence):
    """Formule d'indice Aliphatique(Aal,Val,Leu,Ile) Glycine est ignoré (chaine= 1 Hydrogen)"""
    if not sequence:
        return 0
    
    totalIndex=0
    
    aliphatic=['Ala','Val','Ile','Leu']

    for AA in sequence: 
        if AA.code in aliphatic:
            totalIndex+=AA.ali_weight

    AI=(totalIndex/len(sequence))*100 #not an average but a -> "Global index!"
    return AI

#Dummy test: 
from AA_BD import AA_DB

sequence=[AA_DB["Ala"],AA_DB["Val"],AA_DB["Trp"],AA_DB["Phe"]]
print(f"Stability AI average: {calculate_stability(sequence)}")


def calculate_binding_affinity(sequence):
    """measures how 'sticky' the seqeunce is as in how much it will bind to a target(receptors, other proteins)
       : Sum of Electrostatics (|Charge|) + Dipoles (polarity Sum)"""
    total_charge=sum(abs(AA.charge) for AA in sequence)
    total_polarity=sum(AA.polarity for AA in sequence)

    return total_charge + total_polarity
    #Why SUM -> BINDING = ADDITIVE property !! each charge/polar residue controbutes to the toal affinity!!! 

    
#Normalization des valeurs propriétes des AA: 
#POURQUOI ? -> Normlisation assure que un changement de 10% en Binding Affinity
#A la meme importance qu'un changement de 10% en Stabilité 


def normalize_aa(AA):
    """Retorune un Dictionnaire des valeurs Normalizés des AA(entre 0.0 et 1.0)"""

    normalized_AA={}

    properties=['hydro','polarity','charge','molecular_mass','ali_weight']

    for p in properties:
        val=getattr(AA,p) #Aagin AA.p wouldn't work..
        min=LIMITS[p]['min']
        max=LIMITS[p]['max']

        normalized_AA[p]=(val-min)/(max-min)

    return normalized_AA

#Test
print(f"Normalized AA Test: {normalize_aa(AA_DB["Trp"])}")


#Sliding Window -> Check pour les 3 AA récents ajoutés 
def sliding_window_check(current_sequence):
    """un Check pour les AAs récemment ajoutés: pruning trés rapide
        Retorune True ou false"""

    if len(current_sequence)<3:
        return True
    
    #Steric hindrance
    
    window=current_sequence[-3:] #3 AA récents 

    print(f"BULKY_AAS: {[repr(x) for x in BULKY_AAS]}")
    print(f"WINDOW CODES: {[repr(aa.code) for aa in window]}")


    bulky_count=sum(1 for AA in window if AA.code in BULKY_AAS) #"Generator Expression wow"

    print(f"DEBUG: bulky_count calculated as: {bulky_count}")

    if bulky_count>2: #100% du window = bulky
        return False #Failed the test -> PRUNE
    
    #Static repulsion 
    charge=[aa.charge for aa in window] #list des charges Ex: [-1 +1 +1]..

    

    #Si meme charge -> répulsion -> PRUNE
    if all(c>0 for c in charge):
        return False
    
    if all(c<0 for c in charge):
        return False

    return True


#Test 
sequence2=[AA_DB["Trp"],AA_DB["Phe"],AA_DB["Tyr"]]
print("----------------------------")
print(sliding_window_check(sequence2))

#BRANCH AND BOUND! 
def branch_and_bound_check(current_sequence,user_score):
    current_score=sum()
    #Ghodwa nkamlou el python el kol !
    #mazzle ken el fonction hethi ou 2 fonctions o5rin mta3 
    #Validation scientifique (comparsion avec d'autres peptides kima Insulin)
    #ou wa7da mta3 DFS :) 





