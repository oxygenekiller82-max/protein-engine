package RealStuffs.AA_PFA.repositories;

import java.util.Optional;

import org.springframework.data.jpa.repository.JpaRepository;

import RealStuffs.AA_PFA.model.User;

//no more JpaRepo dang.. NO IT IS BACK wth..
public interface UserRepository extends JpaRepository<User, Long> {
	Optional<User> findByEmail(String username);
	
	Boolean existsByEmail(String username);
	//sign up -> user exists ?
	
	
}
