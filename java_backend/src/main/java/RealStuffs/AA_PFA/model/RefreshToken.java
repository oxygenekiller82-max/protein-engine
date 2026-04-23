package RealStuffs.AA_PFA.model;

import java.time.Instant;

import com.fasterxml.jackson.annotation.JsonIgnore;

import jakarta.persistence.Column;
import jakarta.persistence.Entity;
import jakarta.persistence.GeneratedValue;
import jakarta.persistence.GenerationType;
import jakarta.persistence.Id;
import jakarta.persistence.JoinColumn;
import jakarta.persistence.OneToOne;
import lombok.Data;

@Entity
@Data
public class RefreshToken {
	@Id
    @GeneratedValue(strategy = GenerationType.AUTO)
    private long id;
	
	//relation avec user = 1-1 -> column + 
	@OneToOne
    @JoinColumn(name = "user_id", referencedColumnName = "id")
	@JsonIgnore
    private User user;
	
	@Column(nullable = false, unique = true)
    private String token;
	
	@Column(nullable = false)
    private Instant expiryDate;
}
