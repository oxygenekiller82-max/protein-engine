from AA_BD import AA_DB
from AA_Math import calculate_hydrophobicity, calculate_molecular_mass, \
calculate_stability, calculate_binding_affinity

import matplotlib.pyplot as plt
import numpy as np


#1 letter code -> 3 letter code
DB_mapping={
        'A': 'Ala',
        'R': 'Arg',
        'N': 'Asn',
        'D': 'Asp',
        'C': 'Cys',
        'Q': 'Gln',
        'E': 'Glu',
        'G': 'Gly',
        'H': 'His', 
        'I': 'Ile',
        'L': 'Leu',
        'K': 'Lys',
        'M': 'Met',
        'F': 'Phe',
        'P': 'Pro',
        'S': 'Ser',
        'T': 'Thr',
        'W': 'Trp',
        'Y': 'Tyr',
        'V': 'Val'
    }


#REVERSE MAPPING: 3letter code -> 1 letter code :
#Dictionary comprehension exists...
reverse_mapping ={v: k for k, v in DB_mapping.items()} #v:k ->old value = new key, old key = new value

#iterating through orignal mapping: for k, v in DB_mapping.items() 


# 1 and 3 letter codes -> point to same objects in AA_BD
for one_letter_code, three_letter_code in DB_mapping.items():
    if three_letter_code in AA_DB:
        AA_DB[one_letter_code] = AA_DB[three_letter_code] #POINEUR!
    else: 
        print(f"Fatal:{three_letter_code} not found in AA_DB.")

def getSequenceData(sequence):
    """returns a peptide's data(hydro,mass,stability,binding_affinity)
        sequence must be a string, use 1 LETTER amino acids, not 3"""
        
    three_letter_code_sequence=[]
    for AA in sequence:
        three_letter_code_sequence.append(AA_DB[AA])

    #calcualtions: our helper functions: 
    hydro=calculate_hydrophobicity(three_letter_code_sequence)
    mass=calculate_molecular_mass(three_letter_code_sequence)
    binding_affinity=calculate_binding_affinity(three_letter_code_sequence)
    stability=round(calculate_stability(three_letter_code_sequence,len(three_letter_code_sequence)),2)

    return {
        "hydro": hydro,
        "mass": mass,
        "binding":binding_affinity,
        "stability":stability,
        "length": len(sequence)
        }
    

print("/_/_/_/_//_/_/_/_/_/")
#print(getSequenceData("TDNEPMQKSLCIHYVFGARTTEDNLP"))

print(getSequenceData("MSTVEEDSDTVTVETVNSVTLTQDTEGNLILHCPQNEADEIDSEDSIEPPHKRLCLSSEDDQSIDDSTPCISVVALPLSENDQSFEVTMTATTEVADDEVTEGTVTQIQILQNEQLDEISPLGNEEVSAVSQAWFTTKEDKDSLTNKGHKWKQGMWSKEEIDILMNNIERYLKARGIKDATEIIFEMSKDERKDFYRTIAWGLNRPLFAVYRRVLRMYDDRNHVGKYTPEEIEKLKELRIKHGNDWATIGAALGRSASSVKDRCRLMKDTCNTGKWTEEEEKRLAEVVHELTSTEPGDIVTQGVSWAAVAERVGTRSEKQCRSKWLNYLNWKQSGGTEWTKEDEINLILRIAELDVADENDINWDLLAEGWSSVRSPQWLRSKWWTIKRQIANHKDVSFPVLIKGLKQLHENQKNNPTLLENKSGSGVPNSNTNSSVQHVQIRVARLEDNTAISSSPMAALQIPVQITHVSSADSPATVDSETITLNSGTLQTFEILPSFHLQPTGTPGTYLLQTSSSQGLPLTLTASPTVTLTAAAPASPEQIIVHALSPEHLLNTSDNVTVQCHTPRVIIQTVATEDITSSISQAELTVDSDIQSSDFPEPPDALEADTFPDEIHHPKMTVEPSFNDAHVSKFSDQNSTELMNSVMVRTEEEISDTDLKQEESPSDLASAYVTEGLESPTIEEQVDQTIDDETILIVPSPHGFIQASDVIDTESVLPLTTLTDPILQHHQEESNIIGSSLGSPVSEDSKDVEDLVNCH"))

#The legendary DMTF1 , the one that inspired this project! 
#'hydro': -0.5335526315789474, 'mass': 98145.35, 'binding': 6921.9, 'stability': 83.74, 'length': 760
#i wonder if our program will even run lmao 760 ??
#



#return % of similarity from the arrays 
def calculate_auc_similarity(target_array,res_array): 

    target_auc = np.trapz(np.abs(target_array))
    res_auc = np.trapz(np.abs(res_array))

    if target_auc == 0:
        return 0.0
    
    similarity = (1 - abs(target_auc - res_auc) / target_auc) * 100
    return round(float(similarity), 2)




#Charts ? 
def compare_sequences(target_seq, res_sequence):
    """returns TWO arrays of  for each property, one for the original sequence
        and one for the DFS result sequence, must input both"""
    #properties: hydro for testing 
    target_hydro = [AA_DB[AA].hydro for AA in target_seq]
    res_hydro= [AA_DB[AA].hydro for AA in res_sequence]

    #Stability
    target_stability = [calculate_stability([AA_DB[AA]],1) for AA in target_seq]
    res_stability= [calculate_stability([AA_DB[AA]],1) for AA in res_sequence]

    #Binding
    target_binding=[calculate_binding_affinity([AA_DB[AA]]) for AA in target_seq]
    res_binding=[calculate_binding_affinity([AA_DB[AA]]) for AA in res_sequence]

    #Mass 
    target_mass=[calculate_molecular_mass([AA_DB[AA]]) for AA in target_seq]
    res_mass=[calculate_molecular_mass([AA_DB[AA]]) for AA in res_sequence]

    #AUC SIMILARITY 
    hydro_sim = calculate_auc_similarity(target_hydro,res_hydro)
    stability_sim = calculate_auc_similarity(target_stability, res_stability)
    molecular_mass_sim = calculate_auc_similarity(target_mass, res_mass)
    binding_sim=calculate_auc_similarity(target_binding, res_binding)



    #plt.figure(figsize=(10, 5))

    #plt.plot(indices, target_binding, label='Target Sequence', color='blue', linewidth=2)
    #plt.plot(indices, res_binding, label='DFS Result', color='orange', linestyle='--', linewidth=2)
    
    #plt.title("Hydropathy Profile Comparison")
    #plt.xlabel("Amino Acid Position")
    #plt.ylabel("Hydropathy Score")
    #plt.legend()
    #plt.grid(True, alpha=0.3)
    
    #plt.show()

    #let's just make it return the 8 arrays.
    return {
        "target_hydro": target_hydro, 
        "res_hydro":res_hydro, 
        "target_stability": target_stability,
        "res_stability": res_stability,
        "target_binding": target_binding,
        "res_binding": res_binding,
        "target_mass": target_mass,
        "res_mass": res_mass,

        "hydro_sim":hydro_sim,
        "stability_sim":stability_sim,
        "molecular_mass_sim":molecular_mass_sim,
        "binding_sim":binding_sim,
    }



#compare_sequences("KERSWFVNDYQLAMPTRICHGAKLDV",
#                  "TDNEPMQKSLCIHYVFGARTTEDNLP")
