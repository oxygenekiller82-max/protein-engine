package RealStuffs.AA_PFA.services;

import java.time.LocalDateTime;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.springframework.beans.factory.annotation.Value;
import org.springframework.core.ParameterizedTypeReference;
import org.springframework.security.core.context.ReactiveSecurityContextHolder;
import org.springframework.stereotype.Service;
import org.springframework.web.reactive.function.client.WebClient;

import RealStuffs.AA_PFA.model.ContraintesBiochimiques;
import RealStuffs.AA_PFA.model.Sequence;
import RealStuffs.AA_PFA.repositories.SequenceRepository;
import RealStuffs.AA_PFA.repositories.UserRepository;
import reactor.core.publisher.Mono;
import reactor.core.scheduler.Schedulers;

@Service
public class SequenceService {
	private final WebClient webClient;
    private final SequenceRepository sequenceRepository;
    private final UserRepository userRepository;
    
    @Value("${app.debug.test}")
    private String debugTest;
    
    //INJECTION DE DEPENDANCES! !! 
    //just add these 2 tools ! 
    public SequenceService(WebClient webClient,
    	   SequenceRepository sequenceRepository,
    	   UserRepository userRepository) {
    	
    	this.webClient = webClient;
        this.sequenceRepository = sequenceRepository;
        this.userRepository = userRepository;
    }
    
    //Mono  = async, will wait for the API's reponse in the background
    public Mono<Map<String,Object>> generateAndSave(Map<String,Object> targets,Integer length,boolean bioSwitch){
    	
    	//SECURITY INSTEAD OF USER ID .. omg this method 
    	return ReactiveSecurityContextHolder.getContext()
    			.map(ctx->ctx.getAuthentication().getName())//name is email ofc
    			
    			.flatMap(email->{//->User OBJECT
    				//again the problem of JPA blocking into reactive.. 
    				//will just wrpa into Mono.fromCallable 
    				return Mono.fromCallable(() -> userRepository.findByEmail(email).orElse(null))
    	                       .subscribeOn(Schedulers.boundedElastic());
    				
    				})
    			.flatMap(user->{ //Omg it STAYS! the user object is acceisble !
    				//called "closure" in java..
    				if (user == null) return Mono.error(new RuntimeException("User not found"));
    		
			    	//we have a dict inside a dict in json.. woops
			    	Map<String,Object> requestBody = new HashMap<>();
			    	// the COMPLETE json map, targets, targetLength,BioSwitch
			    	requestBody.put("targets", targets); // TARGETS DICTIONNAIRE!
			    	requestBody.put("target_length", length);
			        requestBody.put("biological_switch", bioSwitch);
			    	
			        
			    	//CALL FLASK API 
			        //WAIT  STEP 1 = generate sequence with null for peoprties just /generate endpoint
			    	return webClient.post()
			    			.uri("/generate") //endpoint 
			    			.bodyValue(requestBody)//convert to json
			    			
			    			.retrieve()//get from python
			    			
			    			//ParameterizedTypeReference.. java FORGETS the types inside the list if the code still going LMAOO
			    			//it fogets the type of the keys
			    			//so tell it .. Map keys = strings, val = objects!
			    			
			    			.bodyToMono(new ParameterizedTypeReference <Map<String,Object>>(){})//CONVERT result to MAP -> keys "function_calls"..
			    			.flatMap(response->{ 
			    				//flatMap = ?? TODO
			    				// 1- SAVE TO DB
			    				
			    				//Sequence seq=new Sequence(); 
			    				
			    				//message de l'API
			    				String message = String.valueOf(response.get("message"));
			    				
			    				//1- sequence list (python has ["Pro","Asp"] .. -> java = "ProAsp"
			    				
			    				@SuppressWarnings("unchecked")
			    				List<String> sequence = (List<String>) response.get("sequence");
			    				
			    				// (sequence != null && !sequence.isEmpty()) {
			    	             //   seq.setPeptideChain(String.join(" ", sequence));
			    	            //} else {
			    	             //   seq.setPeptideChain("API Error: " + message);
			    	            //}
			    				
			    				if (sequence == null || sequence.isEmpty()) {
			    					Sequence errorSeq = new Sequence();
			    					//error mesage save! 
			    					errorSeq.setPeptideChain("API Error: " + message);
			    			        errorSeq.setCreatedAt(LocalDateTime.now());
			    			        errorSeq.setUser(user);
			    			        
			    			        //FROM OROGINAL PARAMTERS! 
			    			        errorSeq.setTargetLength(length);
			    			        errorSeq.setIsBiological(bioSwitch);
			    			        
			    			        sequenceRepository.save(errorSeq);
			    			        return Mono.just(response);
			    			    }
			    				
			    				//STEP 2 -> map to fill in the 4 properties 
			    				Map<String, Object> statsBody = new HashMap<>();
			                    statsBody.put("sequence", sequence);
			                    
			    				//2-STATS -> TO DB AS WELL now
			    				return webClient.post()
			    						.uri("/get_peptide_stats")
			                            .bodyValue(statsBody)
			                            .retrieve()
			                            .bodyToMono(new ParameterizedTypeReference<Map<String,Object>>(){})
			                            //NOO WE HAVE A CONSTRAINTES TABLE !! 
			                            .map(extraStats -> {
			                            	//STEP 3 -> save to DB , with nulls for peroprties
			                            	//parent sequence
			                            	Sequence seq = new Sequence();
			                                seq.setPeptideChain(String.join(" ", sequence));
			                                
			                                //the NON proprties first..
			                                seq.setIsBiological(bioSwitch);
			                                seq.setTargetLength(length);
			                                seq.setCreatedAt(LocalDateTime.now());
			                                seq.setUser(user);
			                                
			                                //contraintes child 
			                                ContraintesBiochimiques results = ContraintesBiochimiques.builder()
			                                		.masseCible(Double.valueOf(String.valueOf(extraStats.get("mass"))))
			                                		.echelleKyteDoolittle(Double.valueOf(String.valueOf(extraStats.get("hydro"))))
			                                		.indiceAliphatique(Double.valueOf(String.valueOf(extraStats.get("stability"))))
			                                		.bindingAffinity(Double.valueOf(String.valueOf(extraStats.get("binding"))))
			                                		// -> seqeunce
			                                		.sequence(seq)
			                                		.build();
			                                
			                                //LINK to 
			                                seq.setContraintes(results);
			                                
			                                //save -> cascade all saves BOTH
			                                sequenceRepository.save(seq);
			                                	
			                                
			                                return response; //HISTORY
			                                //WITH THIS APPROACH
			                                //-> THE HISTTORY WILL LIVE IN THE BROWSER RAM! 
			                                //BUT BUT BUT  to keep it while user opens naothe tab 
			                                //-> must be saved in a "BehaviorSubject"
			                            	
			                            });
			    			});
    			});
		    }
    
  //SECOND API END POINT -> the arrays for charts ! 
	public Mono<Map<String, Object>> getComparisonData(String targetSeq, String generatedSeq){
		Map<String, Object> requestBody = new HashMap<>();
		//keys
		requestBody.put("target_seq", targetSeq);
	    requestBody.put("generated_seq", generatedSeq);
		
	    //no processing (ithink) just put the result a BRIDGE 
	    return webClient.post()
	    		.uri("/generate_arrays")
	    		.bodyValue(requestBody) // -> converts to JSON so python udnetstands it 
	    		//-> PUTS IT IN THE BODY OF THE REQUEST
	    		.retrieve()
	    		.bodyToMono(new ParameterizedTypeReference<Map<String, Object>>() {});
	    		//bodyToMove -> when spring receives from python 
	    		//-> TURNS IT FROM JSON -> JAVA object !
	    
	    		//paramtrizedType again...
	    		//type erasure again.. once java is compiled -> FORGETS THE TYPES INSIDE LIST OR MAP !!
	    		//so it sees this Map<String, Object> AS Map<Unknown, Unknown> !!!
	    		//so response.get("sequence") .. welp won't work 
	    		//ParameterizedTypeReference FORCES java to REMMEBER THE TYPE !!
	}
	
	
	//Get the generated sequence info!!
	public Mono<Map<String,Object>> getSequenceData(List<String> seq){
		Map<String,Object> requestBody = new HashMap<>();
		requestBody.put("sequence",seq);
		System.out.println("PROPS DEBUG: " + debugTest);
		
		//flask 
		return webClient.post()
				.uri("/get_peptide_stats")
				.bodyValue(requestBody)
				.retrieve()
				.bodyToMono(new ParameterizedTypeReference<Map<String,Object>>(){});
		
		
		

	}
	

}


