package RealStuffs.AA_PFA.model;

import java.time.LocalDateTime;

import jakarta.persistence.Column;
import jakarta.persistence.Entity;
import jakarta.persistence.GeneratedValue;
import jakarta.persistence.GenerationType;
import jakarta.persistence.Id;
import jakarta.persistence.JoinColumn;
import jakarta.persistence.ManyToOne;
import jakarta.persistence.PrePersist;
import jakarta.persistence.Table;
import lombok.ToString;
@Entity
@Table(name="sequences")
public class Sequence {
	@Id
    @GeneratedValue(strategy = GenerationType.IDENTITY)
    private Long id;
	
	@Column(columnDefinition = "TEXT")
    private String peptideChain;
	
	private Double mass;
    private Double hydrophobicity;
    private Double stability;
    private Double bindingAffinity;
    
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
	
	public Long getId() {
		return id;
	}

	public void setId(Long id) {
		this.id = id;
	}

	public String getPeptideChain() {
		return peptideChain;
	}

	public void setPeptideChain(String peptideChain) {
		this.peptideChain = peptideChain;
	}

	public Double getMass() {
		return mass;
	}

	public void setMass(Double mass) {
		this.mass = mass;
	}

	public Double getHydrophobicity() {
		return hydrophobicity;
	}

	public void setHydrophobicity(Double hydrophobicity) {
		this.hydrophobicity = hydrophobicity;
	}

	public Double getStability() {
		return stability;
	}

	public void setStability(Double stabilty) {
		this.stability = stabilty;
	}

	public Double getBindingAffinity() {
		return bindingAffinity;
	}

	public void setBindingAffinity(Double bindingAffinity) {
		this.bindingAffinity = bindingAffinity;
	}

	public Integer getFunctionCalls() {
		return functionCalls;
	}

	public void setFunctionCalls(Integer functionCalls) {
		this.functionCalls = functionCalls;
	}

	public Integer getBranchesPruned() {
		return branchesPruned;
	}

	public void setBranchesPruned(Integer branchesPruned) {
		this.branchesPruned = branchesPruned;
	}

	public Boolean getIsBiological() {
		return isBiological;
	}

	public void setIsBiological(Boolean isBiological) {
		this.isBiological = isBiological;
	}

	public LocalDateTime getCreatedAt() {
		return createdAt;
	}

	public void setCreatedAt(LocalDateTime createdAt) {
		this.createdAt = createdAt;
	}

	public User getUser() {
		return user;
	}

	public void setUser(User user) {
		this.user = user;
	}

	
    public Sequence() {
		super();
	}

	@PrePersist
    protected void onCreate() {
        this.createdAt = LocalDateTime.now();
    }

	public Integer getTargetLength() {
		return targetLength;
	}

	public void setTargetLength(Integer targetLength) {
		this.targetLength = targetLength;
	}

}
