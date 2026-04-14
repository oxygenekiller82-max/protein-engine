package RealStuffs.AA_PFA.model;

import jakarta.persistence.Entity;
import jakarta.persistence.GeneratedValue;
import jakarta.persistence.GenerationType;
import jakarta.persistence.Id;
import jakarta.persistence.JoinColumn;
import jakarta.persistence.OneToOne;
import lombok.AllArgsConstructor;
import lombok.Builder;
import lombok.Data;
import lombok.NoArgsConstructor;
import lombok.ToString;

@Entity 
@Data 
@Builder 
@NoArgsConstructor  
@AllArgsConstructor
@ToString(exclude = "sequence")
//bidrectionnele ..so printing one pritnts the other 
public class ContraintesBiochimiques {
	@Id 
	@GeneratedValue(strategy = GenerationType.IDENTITY)
	private Long id; 
	
	private Double masseCible; 
	private Double echelleKyteDoolittle; 
	private Double indiceAliphatique; 
	
	//FORGOT IT ! 
	private Double bindingAffinity; 
	
	//relation avec Sequence -> une sequence a UNE SEULE table de contraintes  1 <-> 1
	
	@OneToOne
	@JoinColumn (name="fk_peptide",referencedColumnName="id")
	//has fk -> it's managing it heh.. 
	private Sequence sequence;
	
}
