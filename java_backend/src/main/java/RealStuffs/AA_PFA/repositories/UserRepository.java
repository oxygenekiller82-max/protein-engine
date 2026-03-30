package RealStuffs.AA_PFA.repositories;

import java.util.Optional;

import org.springframework.data.jpa.repository.JpaRepository;

import RealStuffs.AA_PFA.model.User;

public interface UserRepository extends JpaRepository<User,Long> {
	Optional<User> findByUsername(String username);
	
	boolean existsByUsername(String username);
	//sign up -> user exists ?
}
