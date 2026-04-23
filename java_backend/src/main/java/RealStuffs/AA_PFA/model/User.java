package RealStuffs.AA_PFA.model;

import java.util.Collection;
import java.util.List;

import org.springframework.security.core.GrantedAuthority;
import org.springframework.security.core.authority.SimpleGrantedAuthority;
import org.springframework.security.core.userdetails.UserDetails;

import com.fasterxml.jackson.annotation.JsonIgnore;

import jakarta.persistence.CascadeType;
import jakarta.persistence.Column;
import jakarta.persistence.Entity;
import jakarta.persistence.GeneratedValue;
import jakarta.persistence.GenerationType;
import jakarta.persistence.Id;
import jakarta.persistence.OneToMany;
import jakarta.persistence.Table;
import lombok.AllArgsConstructor;
import lombok.Builder;
import lombok.Data;
import lombok.NoArgsConstructor;
import lombok.ToString;

@Entity 	
@Table(name="users")
@Builder
@Data
@AllArgsConstructor
@NoArgsConstructor
public class User implements UserDetails {
	@Id 
	@GeneratedValue(strategy = GenerationType.IDENTITY)
	private Long id; 
	
	public Long getId() {
		return id;
	}

	public void setId(Long id) {
		this.id = id;
	}

	public String getActualUsername() {
		//oof..
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

	
	public List<Sequence> getSavedSequences() {
		return savedSequences;
	}

	public void setSavedSequences(List<Sequence> savedSequences) {
		this.savedSequences = savedSequences;
	}

	//UserDetails wants to override before working.. okay.. 
	@Override
	//good to get Authorities ngl anwyays 
    public Collection<? extends GrantedAuthority> getAuthorities() {
        // everyone is USER role for now //TODO brag
        return List.of(new SimpleGrantedAuthority("ROLE_USER"));
    }
	
	@Override
	//this for security 
    public String getUsername() {
        return this.email;
    }
	
	
	@Override
    public boolean isAccountNonExpired() { return true; }

    @Override
    public boolean isAccountNonLocked() { return true; }

    @Override
    public boolean isCredentialsNonExpired() { return true; }

    @Override
    public boolean isEnabled() { return true; }
    
	

	@Column(unique = true, nullable = false)
    private String username;
	
	@Column(unique=true,nullable=false)
	private String password; //-> BCrypt pour hashage
	
	@Column(unique = true, nullable = false)
    private String email;
	
	// relation avec seqeunces -> User has many sequences
	
	@OneToMany(mappedBy = "user", cascade = CascadeType.ALL) 
	//cascade all -> pas de sqeunces orphan (elur User supprimée) -> seront supprimés
	@ToString.Exclude //User -> prints sequences -> sequences prints User -> CRASH this prevents it 
	@JsonIgnore 
    private List<Sequence> savedSequences;
	
}
