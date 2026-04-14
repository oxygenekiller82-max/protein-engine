package RealStuffs.AA_PFA.model;

import java.time.LocalDateTime;

import jakarta.persistence.CascadeType;
import jakarta.persistence.Column;
import jakarta.persistence.Entity;
import jakarta.persistence.GeneratedValue;
import jakarta.persistence.GenerationType;
import jakarta.persistence.Id;
import jakarta.persistence.JoinColumn;
import jakarta.persistence.ManyToOne;
import jakarta.persistence.OneToOne;
import jakarta.persistence.PrePersist;
import jakarta.persistence.Table;
import lombok.AllArgsConstructor;
import lombok.Builder;
import lombok.Data;
import lombok.NoArgsConstructor;
import lombok.ToString;


@Entity
@Table(name="sequences")
@Builder 
@Data
@NoArgsConstructor  
@AllArgsConstructor
public class Sequence {
	@Id
    @GeneratedValue(strategy = GenerationType.IDENTITY)
    private Long id;
	
	@Column(columnDefinition = "TEXT")
    private String peptideChain;
	
 
    private Integer functionCalls;
    private Integer branchesPruned;
    private Boolean isBiological; 
    
    private Integer targetLength;

    private LocalDateTime createdAt;
    
    //relation avec user : many sequences same user
    @ManyToOne
    @JoinColumn(name = "user_id")
    @ToString.Exclude
    private User user;
	
	
	@PrePersist
    protected void onCreate() {
        this.createdAt = LocalDateTime.now();
    }

	
	//relation avec Contraintes table
	@OneToOne(mappedBy="sequence",cascade=CascadeType.ALL,orphanRemoval=true)
	//why cascade all ? 
	//-> sequenceRepository.save(seq)
	//-> Sees in the seq object -> Contraintes object attached !! 
	//-> AUTOMATICALLTY saves the sequence to get an ID first 
	//-> saved the cotnraintes table using the fk from that just saved seq object Wow
	private ContraintesBiochimiques contraintes;

}
