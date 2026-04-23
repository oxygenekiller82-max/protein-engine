package RealStuffs.AA_PFA.model;

import lombok.AllArgsConstructor;
import lombok.Builder;
import lombok.Data;
import lombok.NoArgsConstructor;


public class SequenceDTO {
    public Long id;
    public String peptideChain;
    public Integer targetLength;
    public Boolean isBiological;
    public String createdAt;
    public Double masseCible;
    public Double echelleKyteDoolittle;
    public Double indiceAliphatique;
    public Double bindingAffinity;

    public SequenceDTO(Long id, String peptideChain, Integer targetLength, Boolean isBiological,
                       String createdAt, Double masseCible, Double echelleKyteDoolittle,
                       Double indiceAliphatique, Double bindingAffinity) {
        this.id = id;
        this.peptideChain = peptideChain;
        this.targetLength = targetLength;
        this.isBiological = isBiological;
        this.createdAt = createdAt;
        this.masseCible = masseCible;
        this.echelleKyteDoolittle = echelleKyteDoolittle;
        this.indiceAliphatique = indiceAliphatique;
        this.bindingAffinity = bindingAffinity;
    }
}