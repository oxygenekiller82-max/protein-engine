class AminoAcid:
    def __init__(self,name,code,hydro,polarity,charge,is_aliphatic,ali_weight,molecular_mass):
        self.name=name
        self.code=code
        self.hydro=hydro
        self.ploarity=polarity
        self.charge=charge
        self.is_aliphatic=is_aliphatic
        self.ali_weight=ali_weight
        self.molecular_mass=molecular_mass


    def __repr_(self):
        """returns AA's Code"""
        return f"AA({self.code})"
    


#Dictionnaire pour stocker les informations sur AA
AA_DB={
    "A": AminoAcid("Alanine","A",1.8,8.1,0,True,1.0,89.1),
    
    
}