package RealStuffs.AA_PFA.services;

import java.time.LocalDateTime;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.springframework.core.ParameterizedTypeReference;
import org.springframework.stereotype.Service;
import org.springframework.web.reactive.function.client.WebClient;

import RealStuffs.AA_PFA.model.Sequence;
import RealStuffs.AA_PFA.model.User;
import RealStuffs.AA_PFA.repositories.SequenceRepository;
import RealStuffs.AA_PFA.repositories.UserRepository;
import reactor.core.publisher.Mono;

@Service
public class SequenceService {
	private final WebClient webClient;
    private final SequenceRepository sequenceRepository;
    private final UserRepository userRepository;
    
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
    public Mono<Sequence> generateAndSave(Long userId, Map<String,Object> targets,Integer length,boolean bioSwitch){
    	
    	//we have a dict inside a dict in json.. woops
    	Map<String,Object> requestBody = new HashMap<>();
    	// the COMPLETE json map, targets, targetLength,BioSwitch
    	requestBody.put("targets", targets); // TARGETS DICTIONNAIRE!
    	requestBody.put("target_length", length);
        requestBody.put("biological_switch", bioSwitch);
    	
    	
    	//CALL FLASK API 
    	return webClient.post()
    			.uri("/generate") //endpoint 
    			.bodyValue(requestBody)//convert to json
    			
    			.retrieve()//get from python
    			
    			//ParameterizedTypeReference.. java FORGETS the types inside the list if the code still going LMAOO
    			//it fogets the type of the keys
    			//so tell it .. Map keys = strings, val = objects!
    			.bodyToMono(new ParameterizedTypeReference <Map<String,Object>>(){})//CONVERT result to MAP -> keys "function_calls"..
    			.map(response->{ 
    				//PYTHON JSON >- Java entity 
    				Sequence seq=new Sequence(); 
    				
    				//message de l'API
    				String message = String.valueOf(response.get("message"));
    				
    				//1- sequence list (python has ["Pro","Asp"] .. -> java = "ProAsp"
    				
    				@SuppressWarnings("unchecked")
    				List<String> sequence = (List<String>) response.get("sequence");
    				if (sequence != null && !sequence.isEmpty()) {
    	                seq.setPeptideChain(String.join(" ", sequence));
    	            } else {
    	                seq.setPeptideChain("API Error: " + message);
    	            }
    				
    				//2-nested, stats map:
    				@SuppressWarnings("unchecked")
    				Map<String,Object> stats = (Map<String,Object>) response.get("stats");
    				if (stats!=null) {
    					//object -> string -> integer
    					seq.setFunctionCalls(Integer.valueOf(String.valueOf(stats.get("function_calls"))));
    	                seq.setBranchesPruned(Integer.valueOf(String.valueOf(stats.get("branches_pruned"))));
    				}
    				
    				//3-bio switch, target length
    				seq.setIsBiological(bioSwitch);
    	            seq.setTargetLength(length);
    	            seq.setCreatedAt(LocalDateTime.now());
    				
   	
    				//Link to the user
    				//find user
    				User owner = userRepository.findById(userId).orElse(null);
    				//link 
    				seq.setUser(owner);
    				return sequenceRepository.save(seq);
    	
    			});
    			
    }//.map = CALLBACK => runs ONLY AFTER puthon send  data 
}


