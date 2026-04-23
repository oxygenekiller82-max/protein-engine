package RealStuffs.AA_PFA.repositories;

import java.util.Optional;

import org.springframework.data.jpa.repository.JpaRepository;
import org.springframework.data.jpa.repository.Modifying;

import RealStuffs.AA_PFA.model.RefreshToken;
import RealStuffs.AA_PFA.model.User;

public interface RefreshTokenRepository extends JpaRepository<RefreshToken, Long> {
	Optional<RefreshToken> findByToken(String token);
	
	@Modifying
    int deleteByUser(User user);
	
	Optional<RefreshToken> findByUser(User user);
}
