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
print(getSequenceData("TDNEPMQKSLCIHYVFGARTTEDNLP"))

#Charts ? 
def compare_sequences(target_seq, res_sequence):
    #properties: hydro for testing 
    target_hydro = [AA_DB[AA].hydro for AA in target_seq]
    res_hydro= [AA_DB[AA].hydro for AA in res_sequence]

    #Stability
    target_stability = [calculate_stability([AA_DB[AA]],1) for AA in target_seq]
    res_stability= [calculate_stability([AA_DB[AA]],1) for AA in res_sequence]

    #AREA UNDER CURVE ???#-> trapezoidal area Auc = area under curve
    target_auc = np.trapz(target_stability)
    dfs_auc = np.trapz(res_stability)
    similarity = (1 - abs(target_auc - dfs_auc) / target_auc) * 100
    print(f"Structural Bio-Similarity: {similarity:.2f}%")

    #let's do for hydro!
    target_hydro_auc = np.trapz(np.abs(target_hydro))
    dfs_hydro_auc = np.trapz(np.abs(res_hydro))
    hydro_similarity = (1 - abs(target_hydro_auc - dfs_hydro_auc) / target_hydro_auc) * 100
    print(f"Hydropathy Bio-Similarity: {hydro_similarity:.2f}%")


    indices = range(len(target_seq))
    plt.figure(figsize=(10, 5))

    plt.plot(indices, target_stability, label='Target Sequence', color='blue', linewidth=2)
    plt.plot(indices, res_stability, label='DFS Result', color='orange', linestyle='--', linewidth=2)
    
    plt.title("Hydropathy Profile Comparison")
    plt.xlabel("Amino Acid Position")
    plt.ylabel("Hydropathy Score")
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.show()

compare_sequences("KERSWFVNDYQLAMPTRICHGAKLDV",
                  "TDNEPMQKSLCIHYVFGARTTEDNLP")


    