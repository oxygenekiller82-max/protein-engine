package RealStuffs.AA_PFA.controllers;

import java.util.Map;

import org.springframework.web.bind.annotation.CrossOrigin;
import org.springframework.web.bind.annotation.PathVariable;
import org.springframework.web.bind.annotation.PostMapping;
import org.springframework.web.bind.annotation.RequestBody;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RestController;

import RealStuffs.AA_PFA.model.Sequence;
import RealStuffs.AA_PFA.services.SequenceService;
import reactor.core.publisher.Mono;

@RestController
@RequestMapping("/api/PFA")//every URl in controller starts with this here!!
@CrossOrigin(origins = "http://localhost:4200")//angular will talk to this hmm
public class SequenceController {
	private final SequenceService service; 
	
	//injection de dépendances par constrcuteur not autowried -> modern
	public SequenceController(SequenceService service) {
		this.service=service;
	}
	
	
	@PostMapping("/generate/{userId}")
	public Mono<Sequence> generate (@RequestBody Map<String,Object> payload,@PathVariable Long userId){
		
		Integer length = (Integer) payload.get("target_length");
	    boolean isBio = (boolean) payload.get("biological_switch");
	    
		return service.generateAndSave(userId,payload,length, isBio);
	}
	//@RequestBody Map<String,Object> targets -> turns the JSON into a java map
	//@PathVariable Long userId -> gets userId from URL -> to save the sequence the {userId} lol
}
