package RealStuffs.AA_PFA.model;

import java.util.List;

import com.fasterxml.jackson.annotation.JsonIgnore;

import jakarta.persistence.CascadeType;
import jakarta.persistence.Column;
import jakarta.persistence.Entity;
import jakarta.persistence.GeneratedValue;
import jakarta.persistence.GenerationType;
import jakarta.persistence.Id;
import jakarta.persistence.OneToMany;
import jakarta.persistence.Table;
import lombok.ToString;

@Entity 	
@Table(name="users")
public class User {
	@Id 
	@GeneratedValue(strategy = GenerationType.IDENTITY)
	private Long id; 
	
	public Long getId() {
		return id;
	}

	public void setId(Long id) {
		this.id = id;
	}

	public String getUsername() {
		return username;
	}

	public void setUsername(String username) {
		this.username = username;
	}

	public String getPassword() {
		return password;
	}

	public void setPassword(String password) {
		this.password = password;
	}

	@JsonIgnore
	public List<Sequence> getSavedSequences() {
		return savedSequences;
	}

	public void setSavedSequences(List<Sequence> savedSequences) {
		this.savedSequences = savedSequences;
	}

	public User() {
		super();
	}

	@Column(unique = true, nullable = false)
    private String username;
	
	@Column(unique=true,nullable=false)
	private String password; //-> BCrypt pour hashage
	
	// relation avec seqeunces -> User has many sequences
	
	@OneToMany(mappedBy = "user", cascade = CascadeType.ALL) 
	//cascade all -> pas de sqeunces orphan (elur User supprimée) -> seront supprimés
	@ToString.Exclude //User -> prints sequences -> sequences prints User -> CRASH this prevents it 
    private List<Sequence> savedSequences;
	
}
