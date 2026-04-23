package RealStuffs.AA_PFA.controllers;

import java.util.Map;

import org.springframework.http.HttpStatus;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.CrossOrigin;
import org.springframework.web.bind.annotation.PostMapping;
import org.springframework.web.bind.annotation.RequestBody;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RestController;

import RealStuffs.AA_PFA.dtos.LoginRequest;
import RealStuffs.AA_PFA.dtos.LoginResponse;
import RealStuffs.AA_PFA.dtos.RegisterRequest;
import RealStuffs.AA_PFA.services.AuthenticationService;
import RealStuffs.AA_PFA.services.JwtService;
import lombok.RequiredArgsConstructor;
import reactor.core.publisher.Mono;

@RestController
@RequestMapping("/api/auth")
@RequiredArgsConstructor
@CrossOrigin(origins = "http://localhost:4200")
public class AuthenticationController {
	//dependence
	private final AuthenticationService service;
	private final JwtService jwtService;
	
	//POST: /register
	@PostMapping("/register")
	public Mono<ResponseEntity<Map<String, String>>> register(@RequestBody RegisterRequest request) {
        return service.Register(request)
        .map(user->ResponseEntity.ok(Map.of(
                "message", "Account created.",
                "status", "success"
            )))
        	.onErrorResume(e -> Mono.just(ResponseEntity
        	            .status(HttpStatus.CONFLICT)
        	            .body(Map.of("message", e.getMessage()))));
    }
	//login jwt 
	
	@PostMapping("/login")
	public Mono<ResponseEntity<LoginResponse>> login(@RequestBody LoginRequest request){
		return service.authenticate(request)
				.map(response->ResponseEntity.ok(response))
				.defaultIfEmpty(ResponseEntity.status(HttpStatus.UNAUTHORIZED).build());
	}
	
	//refreshing the token: 
	@PostMapping("/refresh")
	public ResponseEntity<?> refresh(@RequestBody Map<String, String> request){
		String refreshToken = request.get( "refreshToken");
		
		try {
	        String newAccessToken = jwtService.refreshAccessToken(refreshToken);
	        
	        return ResponseEntity.ok(Map.of(
	                "accessToken", newAccessToken,
	                "refreshToken", refreshToken
	            ));
		}catch(RuntimeException e) {
			//why runtime exception tho..
			//->It bubbles up the exceptions! 
			// + s^ring is designed for runtime exceptions :) 
			//automatically turns then inot 500 error 
			return ResponseEntity.status(HttpStatus.UNAUTHORIZED).body(e.getMessage());
			
		}
	}
	
}
