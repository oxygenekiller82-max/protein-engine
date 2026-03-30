package RealStuffs.AA_PFA.repositories;


import java.util.List;

import org.springframework.data.jpa.repository.JpaRepository;

import RealStuffs.AA_PFA.model.Sequence;

public interface SequenceRepository extends JpaRepository<Sequence, Long> {
	List<Sequence> findByUserId(Long userId); //get sequences 
	
	long countByUserId(Long userId);// how many sequences 
	
	List<Sequence> findByUserAndIsBiological (Long userId ,boolean isBiological); 
	//how many bioligical sequences // JPA sees the AND wow
	
	void deleteByIdAndUserId(Long id, Long userId); // maybe needed
}
