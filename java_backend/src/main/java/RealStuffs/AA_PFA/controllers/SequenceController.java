package RealStuffs.AA_PFA.controllers;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.springframework.security.core.context.ReactiveSecurityContextHolder;
import org.springframework.web.bind.annotation.CrossOrigin;
import org.springframework.web.bind.annotation.GetMapping;
import org.springframework.web.bind.annotation.PostMapping;
import org.springframework.web.bind.annotation.RequestBody;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RestController;

import RealStuffs.AA_PFA.model.Sequence;
import RealStuffs.AA_PFA.model.SequenceDTO;
import RealStuffs.AA_PFA.model.User;
import RealStuffs.AA_PFA.repositories.SequenceRepository;
import RealStuffs.AA_PFA.repositories.UserRepository;
import RealStuffs.AA_PFA.services.SequenceService;
import lombok.RequiredArgsConstructor;
import reactor.core.publisher.Flux;
import reactor.core.publisher.Mono;
import reactor.core.scheduler.Schedulers;

@RestController
@RequestMapping("/api/PFA")//every URl in controller starts with this here!!
@CrossOrigin(origins = "http://localhost:4200")//angular will talk to this hmm
@RequiredArgsConstructor
public class SequenceController {
	private final SequenceService service; 
	private final UserRepository userRepo;
	private final SequenceRepository sequenceRepo;
	
	
	//injection de dépendances par constrcuteur not autowried -> modern
	//public SequenceController(SequenceService service) {
	//	this.service=service;
	//}
	
	//SECURITY context and UserDetails gives who's connected lol 
	//NO need for the  /{id} path variable that's really neat! 
	
	@PostMapping("/generate")
	public Mono<Map<String, Object>> generate (@RequestBody Map<String,Object> payload){
		@SuppressWarnings("unchecked")
		
		//INNER STUFF IN THE JSON
		Map<String, Object> targetsMap = (Map<String, Object>) payload.get("targets");
		
		//OUTER STUFF
		Integer length = (Integer) payload.get("target_length");
	    boolean isBio = (boolean) payload.get("biological_switch");
	    
		return service.generateAndSave(targetsMap,length, isBio);
	}
	//@RequestBody Map<String,Object> targets -> turns the JSON into a java map
	//@PathVariable Long userId -> gets userId from URL -> to save the sequence the {userId} lol
	
	@PostMapping("/compare")
	public Mono<Map<String,Object>> compareSequences(@RequestBody Map<String, String> payload){
		String target=payload.get("target_seq");
		String generated = payload.get("generated_seq");
		//sevice function ->  Mono<Map<String, Object>> getComparisonData
		return service.getComparisonData(target, generated);
		
	}
	
	@GetMapping("/my-sequences")
	//DTO :( otherse won't work 
	
	public Flux<SequenceDTO> getUserSequences() {
	    
	    // 1. Reach into the Reactive Security Context
	    return ReactiveSecurityContextHolder.getContext()
	        .map(securityContext -> securityContext.getAuthentication().getName()) // gets the email from the token!
	        .flatMapMany(email -> 
	        
	            // BLOCKING JPA -> wrappied into mono from callable 
	        
	            Mono.fromCallable(() -> {
	                User user = userRepo.findByEmail(email)
	                        .orElseThrow(() -> new RuntimeException("User not found"));
	                
	                List<Sequence> sequences = sequenceRepo.findByUser(user);
	                System.out.println("DEBUG user: " + user.getId() + " | sequences found: " + sequences.size()); 
	                
	                return sequences.stream()
                    .map(s -> new SequenceDTO(
                        s.getId(),
                        s.getPeptideChain(),
                        s.getTargetLength(),
                        s.getIsBiological(),
                        s.getCreatedAt() != null ? s.getCreatedAt().toString() : null,
                        s.getContraintes() != null ? s.getContraintes().getMasseCible() : null,
                        s.getContraintes() != null ? s.getContraintes().getEchelleKyteDoolittle() : null,
                        s.getContraintes() != null ? s.getContraintes().getIndiceAliphatique() : null,
                        s.getContraintes() != null ? s.getContraintes().getBindingAffinity() : null
                    ))
                    .collect(Collectors.toList());
	                
	            }).subscribeOn(Schedulers.boundedElastic()) 
	            .flatMapMany(Flux::fromIterable));
	            //.doOnError(e -> System.err.println("ERROR IN MY-SEQUENCES: " + e.getClass().getName() + " -> " + e.getMessage()))
	            //.doOnSuccess(result -> System.out.println("DEBUG result size: " + (result != null ? result.size() : "NULL"))));
	        
	}
	@PostMapping("/sequence_stats")
	public Mono<Map<String, Object>> getSequenceStats(@RequestBody Map<String, List<String>> payload) {
	    List<String> sequence = payload.get("sequence");
	    return service.getSequenceData(sequence);
	}
	
	
	
}
