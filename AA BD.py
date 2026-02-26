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
    "Ala": AminoAcid("Alanine","Ala",1.8,8.1,0,True,1.0,89.1),
    "Gly": AminoAcid("Glycine","Gly",-0.4,9.0,0,True,0,75.07),
    "Ser": AminoAcid("Serine","Ser",-0.8,9.2,0,False,1.0,105.09),
    "Pro": AminoAcid("Proline","Pro",-1.6,8.0,0,False,0,115.13),
    "Val": AminoAcid("Valine","Val",4.2,5.9,0,True,2.9,117.15),
    "Thr": AminoAcid("Threonine","Thr",-0.7,8.6,0,False,0,119.12),
    "Cys": AminoAcid("Cysteine","Cys",2.5,5.5,0,False,0,121.16),
    "Leu": AminoAcid("Leucine","Leu",3.8,4.9,0,True,3.9,131.18),
    "Ile": AminoAcid("Isoleucine","Ile",4.5,5.2,0,True,3.9,131.18),
    "Asn": AminoAcid("Asparagine","Asn",-3.5,11.6,0,False,0,132.12),
    "Asp": AminoAcid("Aspartic acid","Asp",-3.5,13.0,-1,False,0,133.11),
    "Gln": AminoAcid("Glutamine","Gln",-3.5,10.5,0,False,0,146.15),
    "Lys": AminoAcid("Lysine","Lys",-3.9,11.3,1,False,0,146.19),
    "Glu": AminoAcid("Glutamic acid","Glu",-3.5,12.3,-1,False,0,147.13),
    "Met": AminoAcid("Methionine","Met",1.9,5.7,0,False,0,149.21),
    "His": AminoAcid("Histidine","His",-3.9,10.4,1,False,0,155.16),
    "Phe": AminoAcid("Phenylalanine","Phe",2.8,5.2,0,False,0,165.19),
    "Arg": AminoAcod("Arginine","Arg",-4.5,10.5,1,False,0,174.20),
    "Tyr": AminoAcid("Tyrosine","Tyr",-1.3,6.2,0,False,0,181.19),
    "Trp": AminoAcid("Tryptophan","Trp",-0.9,5.4,0,False,0,204.23),
    
    
}