package RealStuffs.AA_PFA.services;



import org.springframework.http.HttpStatus;
import org.springframework.security.crypto.bcrypt.BCryptPasswordEncoder;
import org.springframework.stereotype.Service;
import org.springframework.web.server.ResponseStatusException;

import RealStuffs.AA_PFA.dtos.LoginRequest;
import RealStuffs.AA_PFA.dtos.LoginResponse;
import RealStuffs.AA_PFA.dtos.RegisterRequest;
import RealStuffs.AA_PFA.model.User;
import RealStuffs.AA_PFA.repositories.UserRepository;
import lombok.RequiredArgsConstructor;
import reactor.core.publisher.Mono;
import reactor.core.scheduler.Schedulers;

@Service
@RequiredArgsConstructor
public class AuthenticationService {
	private final JwtService jwtService;
	
	private final UserRepository userRepo;
	private final BCryptPasswordEncoder passwordEncoder;
	//private final UserMapper mapper;
	//user mapper...  but Angualar send JSON..
	
	public Mono<User> Register(RegisterRequest request) {
		//exists ?
			//cannot do .isPresent in Reactive.. wth 
		return Mono.fromCallable(() -> userRepo.existsByEmail(request.getEmail()))
				.flatMap(exists->{
					
					if(exists) {
						//return Mono.error(new RuntimeException("This Email cannot be used for registering."));
						//TODO how to send these to angular.. hmm
						//LIKE THIS Mono.error
						return Mono.error(new ResponseStatusException(HttpStatus.CONFLICT, "Email already exists"));
					}
					
					//if not just proceed this is reactive 
					String generatedUsername = request.getFirstName() + " " + request.getLastName();

					User newUser = User.builder()
			                .email(request.getEmail())
			                .username(generatedUsername)
			                .password(passwordEncoder.encode(request.getPassword()))
			                .build();
					
					return Mono.fromCallable(() -> userRepo.save(newUser));
				})
				.subscribeOn(Schedulers.boundedElastic()); //WHAT THIS IS VERY IMPORTANTN APARNETLY
		}
		
		
	
		
		
	

	
	//and the Login JWT reactive this time joruney begins.. 
	public Mono<LoginResponse> authenticate(LoginRequest request) {
		//email -> password hash matches -> else excpetion
		
		return Mono.justOrEmpty(userRepo.findByEmail(request.getEmail()))
			
			.filter(user -> passwordEncoder.matches(request.getPassword(), user.getPassword()))
			
	        .map(user -> {
	        	String access = jwtService.generateToken(user); // = 15 mins
	        	String refresh = jwtService.generateRefreshToken(user); // will do 7 days ig
	        	
	            return LoginResponse.builder()
	            		.accessToken(access)
	            		.refreshToken(refresh)
	            		.username(user.getActualUsername())
	            		.build();
	            		
        })
        .switchIfEmpty(Mono.error(new ResponseStatusException(HttpStatus.UNAUTHORIZED)));
	}

}
