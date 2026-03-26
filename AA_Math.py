#Formulas file 
import math 
from AA_BD import LIMITS
from AA_BD import BULKY_AAS
import random
import csv

def calculate_hydrophobicity(sequence):
    """ 1/n * Sum(Hi): returns hydrophobicity result of the list of AA"""
    if len(sequence)==0: return 0 
    total_hydro=sum(aa.hydro for aa in sequence)
    return total_hydro/len(sequence)

#Why not, molecular mass slider ?

def calculate_molecular_mass(sequence):
    """Sum of all AA molecular mass"""
    return sum(AA.molecular_mass for AA in sequence)


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
def calculate_stability(sequence,user_target_length):
    """Formule d'indice Aliphatique(Aal,Val,Leu,Ile) Glycine est ignoré (chaine= 1 Hydrogen)
        UPDATED: divide by the OVERALL length not just the current_sequence length
        otherwise during DFS the average will be based on the current_sequence of AAs only ! 
        which will be a very high not accurate average!!!!!!!!
    """
    if not sequence:
        return 0
    
    totalIndex=0
    
    aliphatic=['Ala','Val','Ile','Leu']

    for AA in sequence: 
        if AA.code in aliphatic:
            totalIndex+=AA.ali_weight

    AI=(totalIndex/user_target_length)*100 #not an average but a -> "Global index!"
    #Simplified version of the X(Ala)+a.X(val).. 
    return AI

#Dummy test: 
from AA_BD import AA_DB

sequence=[AA_DB["Ala"],AA_DB["Val"],AA_DB["Trp"],AA_DB["Phe"]]
print(f"Stability AI average: {calculate_stability(sequence,len(sequence))}")


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

    bulky_count=sum(1 for AA in window if AA.code in BULKY_AAS) #"Generator Expression wow"

    #print(f"bulky_count: {bulky_count}")

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
#binding affinity = fonction de : Polarité + charge -> Max score pour une séquence de longeur K= 
#AA qui a max polarité + charge! -> à déterminer: 

def best_AA_binding_affinity():
    max_score=0
    last_AA_checked=None

    for AA in AA_DB.values():
        #1-Normaliser
        normalized_AA = normalize_aa(AA)
        #2-calculer average
        current_score=(normalized_AA['polarity']+abs(AA.charge))/2
        #3-compare
        if current_score>max_score:
            max_score=current_score
            last_AA_checked=AA

    #print(f"MAX score stabilité (polarité + charge) {max_score}")
    #print(f"AA qui contribue le MAX au stabilité:{last_AA_checked}")
    #FOR NOW with our sources -> it's ... Aspartic Acid! ! 'Asp' ! takes the crown! 
      
    return last_AA_checked

#Et bine sur pour le min bound -> AA qui contribue le moints pour binding affinity
def worst_AA_binding_affinity():
    min_score=float('inf')
    last_AA_checked=None

    for AA in AA_DB.values():
        #1-Normaliser
        normalized_AA = normalize_aa(AA)
        #2-calculer average
        current_score=(normalized_AA['polarity']+abs(AA.charge))/2
        #3-compare
        if current_score<min_score:
            min_score=current_score
            last_AA_checked=AA

    #print(f"MIN score stabilité (polarité + charge) {min_score}")
    #print(f"AA qui contribue le MIN au stabilité:{last_AA_checked}")
    #FOR NOW with our sources -> it's ... Aspartic Acid! ! 'Asp' ! takes the crown! 
      
    return last_AA_checked
    
print(best_AA_binding_affinity())
print(worst_AA_binding_affinity())


def get_branch_and_bound(current_sequence,target_length):
    """Input: -current_sequence: current AA list,
              -user_score for EACH property: a Dict, 
              -user sequence target length: Integer.
        Output: Dict of all max and min scores(val) for EACH property(key)
    """
    properties=['hydro','molecular_mass','stability','binding_affinity']

    #K=remaining length for the sequence = positions libres de AA
    k=target_length-len(current_sequence)

    branch_bounds={} #Result: like this; {'hydro':{MAX_theoretical_score,MIN_theoretical_score}}

    #Partie 1: Hydro + masse molécualaire -> K*best_score 
    branch_bounds['hydro']={
        'min': calculate_hydrophobicity(current_sequence)+(k*LIMITS['hydro']['min'])/target_length,
        'max': calculate_hydrophobicity(current_sequence)+(k*LIMITS['hydro']['max'])/target_length
    }

    #branch_bounds['molecular_mass']={
        #'min': calculate_molecular_mass(current_sequence)+(k*LIMITS['molecular_mass']['min']),
        #'max': calculate_molecular_mass(current_sequence)+(k*LIMITS['molecular_mass']['max']),
    #}

    #Unrealistic so.. 
    realistic_min_avg = 100.0  
    realistic_max_avg = 150.0  

    branch_bounds['molecular_mass'] = {
        'min': calculate_molecular_mass(current_sequence) + (k * realistic_min_avg),
        'max': calculate_molecular_mass(current_sequence) + (k * realistic_max_avg),
    }

    #Partie2: Stability + binding_affinity... -> stability(curreny_sequqence) + MAAX Aliphatic Weight index for the remaining K positions

    branch_bounds['stability']={
        'min': calculate_stability(current_sequence,target_length)+(k*LIMITS['ali_weight']['min']/target_length)*100,
        'max': calculate_stability(current_sequence,target_length)+(k*LIMITS['ali_weight']['max']/target_length)*100 
    }

    
    branch_bounds['binding_affinity']={
        'min': calculate_binding_affinity(current_sequence) + k*calculate_binding_affinity([worst_AA_binding_affinity()]),
        'max': calculate_binding_affinity(current_sequence) + k*calculate_binding_affinity([best_AA_binding_affinity()])
    }

    return branch_bounds

#Testing branch and bound because THIS TOOK AS MUCH TIME AS ALL OTHER FUNCTIONS COMBINED ???
branch_bound_test_sequence=[
    AA_DB['Asp'], AA_DB['Leu'], AA_DB['Asp'], 
    AA_DB['Ala'], AA_DB['Glu'], AA_DB['Leu'], 
    AA_DB['Arg'], AA_DB['Gly'], AA_DB['Lys'], AA_DB['Leu']
]
print("---- BRANCH --- AND --- BOUND YEEHAWW")
print(get_branch_and_bound(branch_bound_test_sequence,51))
print("-------------------------------------")
#WOW


def is_sequence_good(current_sequence,branch_bounds,user_targets):
    """a big True/False check: tells DFS Continue or PRUNE!"""

    #SLIDING WINDOW! 
    if not sliding_window_check(current_sequence):
        return False

    #masse moléculaire: masse lightest already heavier than user target -> PRUNE
    #or masse heaviest already lighter than user target -> PRUNE
    if branch_bounds['molecular_mass']['min']>user_targets['mass_max'] or  branch_bounds['molecular_mass']['max'] < user_targets['mass_min']: 
        return False
    
    #Hydro: 
    if branch_bounds['hydro']['max'] < user_targets['hydro_min'] or branch_bounds['hydro']['min'] > user_targets['hydro_max']:
        return False
    
    #binding affinity: max affinity possible < min de user -> PRUNE
    if branch_bounds['binding_affinity']['max'] < user_targets['binding_min']:
        return False
    
    #Stability : max stabilié possible < min de user -> PRUNE
    if branch_bounds['stability']['max'] < user_targets['stability_min']:
        return False
    
    return True #PASSES!


stats_for_nerds={"branches Pruned":0,
                 "function calls": 0
                }
DFS_history=[]

def add_to_history(DFS_history,action,current_sequence,AA_added="None"):
    """capture de l'algorithme, l'ajout aréps Ajout/Prune action"""
    frame={
        "step":stats_for_nerds["function calls"],
        "action":action,
        "sequence":current_sequence.copy(), #LISTES par addresse.. they will chagen ->.copy = snapshot won't change omg
        "last_AA_added":AA_added,
        "level":len(current_sequence) 
    }
    DFS_history.append(frame)


#GREEEDYYYY 
def user_target_weights(user_targets,target_length):
    """calcualtes the most demanding user targets PER AA -> those the closest to the max possible in the AA_DB 
    = the strictest user_targets"""
    #PER AA the MAX: 
    max_binding_per_AA = (LIMITS['polarity']['max'])+ abs(LIMITS['charge']['max']) # should be 14
    max_stability_per_AA=(LIMITS['ali_weight']['max'])*100


    binding_pressure=(user_targets['binding_min']/target_length) / max_binding_per_AA
    stability_pressure=user_targets['stability_min']/max_stability_per_AA

    user_hydro_range=user_targets['hydro_max']-user_targets['hydro_min']
    hydro_range=LIMITS['hydro']['max'] - LIMITS['hydro']['min']
    hydro_pressure = hydro_range / max(0.1, user_hydro_range) #DIVISON BY 0 DANG ITTT

    user_mass_range=user_targets['mass_max']-user_targets['mass_min']
    mass_range=LIMITS['molecular_mass']['max'] - LIMITS['molecular_mass']['min']
    mass_pressure = mass_range / max(0.1, user_mass_range) 


    return {
        'binding': max(0.1, binding_pressure),
        'stability': max(0.1, stability_pressure),
        'hydro': max(0.1, hydro_pressure),
        'molecular_mass': max(0.1, mass_pressure)
    }


def score_AA(AA,user_targets,target_length,user_target_weights): 
    """gives a single normlaized score for an AA based on how well it fits the user targets"""
    #all calculations 
    #contribution of this AA 
    stability_score=calculate_stability([AA],target_length)
    binding_score=calculate_binding_affinity([AA])

    ideal_hydro_score=(user_targets['hydro_min']+user_targets['hydro_max'])/2 #Average
    hydro_dist = abs(AA.hydro-ideal_hydro_score) #distance from the ideal from our AA score
    #same for mass -> but per AA!!!  BOTH MUST MINIMIZE NOT MASMIZE O:

    ideal_mass_per_aa = (user_targets['mass_min'] + user_targets['mass_max']) / (2 * target_length)
    mass_dist = abs(AA.molecular_mass - ideal_mass_per_aa)

    #FINALLY THE SCORE UGH
    score=((stability_score*user_target_weights['stability'])+ \
           (binding_score*user_target_weights['binding']) - \
           (hydro_dist*user_target_weights['hydro']) - \
           (mass_dist * user_target_weights['molecular_mass'])
           )
    
    return score



  
def best_AA(user_targets,target_length,user_target_weights,current_sequence):
    """returns the BEST AAs (sorted AA_DB) to add in the sequence DFS call : Greedy algorithm
        added: penalty for choooisng the same one over and voer again """
    DB=list(AA_DB.values())

    #History 
    history=[aa.code for aa in current_sequence[:]]

   # original Score - (abs( Score) * Penalty multiplier)
    DB.sort(key=lambda aa: score_AA(aa,user_targets,target_length,user_target_weights)-
            (abs(score_AA(aa, user_targets, target_length, user_target_weights)) *(2.0** history.count(aa.code) -1.0)),
            reverse=True)
    #50 * number of occurences subtracted: the penalty!

    return DB


#FINAL RESULT:
found_AA=[]
sample_targets={
    'hydro_min': -40.0,
    'hydro_max': 0.0,
    'mass_min': 3400,
    'mass_max': 3600,
    'stability_min': 92.0,
    'binding_min': 90.0
}
target_length=30
#


weights = user_target_weights(sample_targets, 51) #!!! must precalcuualte before the DFS !!!!

def DFS(current_sequence,target_length,user_targets,res):
    """Input: user_targets: a dict of user preferances for the generated sequence(min_max of Hydro/mass and min stability/binding_affinity,
              target_length: Integer of desired sequence length,
              current_sequence should be an empty list [] to build upon
              res=empty list to store the result
        Output= the list res
    """

    #DFS = réscursivité:
    #base case: déja seqeunce de logneur target_len trouvée


    stats_for_nerds["function calls"]+=1
    print(f"NODES VISITIED: {stats_for_nerds["function calls"]}")

    if(len(res))>0: #is res the result final list done -> stop
        return True

    if len(current_sequence)==target_length:
        res.append(list(current_sequence)) #ALL VALID LISTS! 
        #Why ? -> res = passage par addresse en python !!
        #when it hits level 51 in the tree -> finds AA -> adds it and good 
        #goes back to lelve 50 -> resumes when it stopped in the loop -> tries the next AA, if 
        # also good -> a new list is formed! 

        clean_seq = [aa.code for aa in current_sequence]

        if len(current_sequence)>0: #current_sequence[-1]) in an empty list ...
            add_to_history(DFS_history,"DONE",clean_seq,clean_seq[-1]) 
        else:
            add_to_history(DFS_history,"DONE",clean_seq,[]) 
            
        return True
    
    branch_bounds=get_branch_and_bound(current_sequence,target_length)
    if not is_sequence_good(current_sequence,branch_bounds,user_targets):
        stats_for_nerds["branches Pruned"]+=1

        clean_seq = [aa.code for aa in current_sequence]

        if len(current_sequence)>0: #current_sequence[-1]) in an empty list ...
            add_to_history(DFS_history,"PRUNE",clean_seq,clean_seq[-1]) 
        else:
            add_to_history(DFS_history,"PRUNE",clean_seq,[]) 
        return False #PRUNE THIS BRANCH
    
    #exploration  NO LONGER RANDOM !!!!
    #randomized_DB=list(AA_DB.values())
    #random.shuffle(randomized_DB) #CHANGE ORDER OF AAs with EVERY STEP ! 
    #solves the it's all Ala problem lol 
    
    for AA in best_AA(user_targets,target_length,weights,current_sequence):
        #ADD AA -> ADD TO HISTORY! 
        
        clean_seq=[aa.code for aa in current_sequence]+[AA.code] #instead of current_sequence+[AA.code]
        add_to_history(DFS_history,"ADD",clean_seq,AA.code)
        if DFS(current_sequence+[AA],target_length,user_targets,res): #ALSO ONLY RETURNS ONE SEQUENCE, no more all sorry
            return True #stops if finds sequence! 
        

        #why a FOR loop ? well.. -> if the NEXT function call fails (returns None which is return btw)
        #-> it will RESUME WHERE IT LEFT OFF IN THE FOR LOOP meaning it will just pick something else and not jsut the program! 

    return False #branch is bad -> STOPS tried all 20AAs (the sorted ones!) but none lead to a valid sequence -> loop just stops
                 #parent of this function call sees it -> now moves to the next AA yaay!!!




def validate_user_targets(user_targets,target_length):
    #Min > max -> DFS won't stop
    if (user_targets['hydro_min']> user_targets['hydro_max']):
        print("Cannot start search: Minimum targets cannot be greater than maximum ones")
        return False
    
    if(user_targets['mass_min']> user_targets['mass_max']):
        print("Cannot start search: Minimum targets cannot be greater than maximum ones")
        return False

    if (target_length<=0):
        print("Cannot start search: length must be strictly greater than zero")
        return False

    #Check -> is the sequence physically possible given the consttraints ?
    #DFS DOES check these but it won't stop if so just moves to the next branch so..

    if user_targets['hydro_max']<(target_length*LIMITS['hydro']['min']):
        print(f"Cannot start search: Hydrophobicity should be at least{target_length*LIMITS['hydro']['min']}")
        return False
    
    if user_targets['hydro_min'] > (target_length * LIMITS['hydro']['max']):
        print(f"Cannot start search: Hydrophobicity for this length is at most {target_length * LIMITS['hydro']['max']}")
        return False
    
    if user_targets['mass_max'] < (target_length * LIMITS['molecular_mass']['min']):
        print(f"Cannot start search: Molecular Mass should be at least{target_length*LIMITS['molecular_mass']['min']}")
        return False
    
    if user_targets['mass_min'] > (target_length * LIMITS['molecular_mass']['max']):
        print(f"Cannot start search: Molecular Mass is for this length is at most {target_length * LIMITS['molecular_mass']['max']}")
        return False
    
    if user_targets['binding_min'] > (target_length*calculate_binding_affinity([best_AA_binding_affinity()])):
        print(f"Cannot start search: Binding Affinity at most is {target_length*calculate_binding_affinity([best_AA_binding_affinity()])}")
        return False
    
    #(k*LIMITS['ali_weight']['min']/target_length)*100, k = target length here..
    if user_targets['stability_min']>(LIMITS['ali_weight']['max']*100):
        print(f"Cannot start search: Stability at most is {LIMITS['ali_weight']['max']*100}")
        return False

    return True


#Test

if (validate_user_targets(sample_targets,target_length)):
    DFS([],target_length,sample_targets,found_AA)

print(f"WEIGHTS:{user_target_weights(sample_targets,target_length)}")

print("****************************************")

if found_AA:
    print(f"AAND WE GOT:{found_AA[0]}")

print(stats_for_nerds)

#Exporting the DFS History to a CSV!!! file: 

def save_DFS_history_to_csv(history,filename="dfs_results.cvs"):
    #COLUMNS HEADERS! 
    fields=["Step","Action","Level","Last_AA","Sequence"]

    with open(filename,mode="w",newline="",encoding="utf-8") as f:
        #HEADER:
        header=header = f"{'STEP':<8} | {'ACTION':<10} | {'LVL':<6} | {'LAST_AA':<10} | {'SEQUENCE'}\n"
        f.write(header)
        f.write("-" * 100 + "\n")

        writer = csv.DictWriter(f,fieldnames=fields)

        writer.writeheader()

        for frame in history: 
            seq_string="-".join([str(aa) for aa in frame['sequence']])
            
            line=(
                f"{frame['step']:<8} | "
                    f"{frame['action']:<10} | "
                    f"{frame['level']:<6} | "
                    f"{str(frame['last_AA_added']):<10} | "
                    f"{seq_string}\n"
                )
            f.write(line)

save_DFS_history_to_csv(DFS_history) 

#user_target_weights(user_targets,target_length):

score_test=score_AA(AA_DB['Ile'],sample_targets,51,user_target_weights(sample_targets,51))
print(f"AA SCORE Ile:{score_test}")

print(score_AA(AA_DB['Trp'],sample_targets,51,user_target_weights(sample_targets,51)))
print(score_AA(AA_DB['Ala'],sample_targets,51,user_target_weights(sample_targets,51)))


# W A I T this is actually insane! 
#best ones so far
#1- 30 mer Defensin 
#2- 60mer Silk Style 
#3- 35mer Red Mamba
#4- 46mer Crambin
#5 51mer inslin 



#Validation scientifique: 



   





