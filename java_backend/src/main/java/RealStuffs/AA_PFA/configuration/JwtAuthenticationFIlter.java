package RealStuffs.AA_PFA.configuration;

import java.util.ArrayList;

import org.springframework.http.HttpHeaders;
import org.springframework.security.authentication.UsernamePasswordAuthenticationToken;
import org.springframework.security.core.context.ReactiveSecurityContextHolder;
import org.springframework.security.core.context.SecurityContext;
import org.springframework.stereotype.Component;
import org.springframework.web.server.ServerWebExchange;
import org.springframework.web.server.WebFilter;
import org.springframework.web.server.WebFilterChain;

import RealStuffs.AA_PFA.services.JwtService;
import lombok.RequiredArgsConstructor;
import reactor.core.publisher.Mono;

@Component
@RequiredArgsConstructor
public class JwtAuthenticationFIlter implements WebFilter {
	
	private final JwtService jwtService;

	@Override
		public Mono<Void> filter(ServerWebExchange exchange, WebFilterChain chain) {
			//getting the header
			String authHeader = exchange.getRequest().getHeaders().getFirst(HttpHeaders.AUTHORIZATION);
			
			//and checking the header 
			if (authHeader == null || !authHeader.startsWith("Bearer ")) {
	            return chain.filter(exchange);
	        }
			//MUST be a value + STARTS WITH BEARER
			
			String jwt = authHeader.substring(7);
			
			//TRY CATCH otherwise bad mlafomred token -> app crash NOO!
			try {
		        String userEmail = jwtService.extractUsername(jwt);
		        
		        return ReactiveSecurityContextHolder.getContext()
		        		.map(SecurityContext::getAuthentication)
		        		.flatMap(auth -> chain.filter(exchange))//authenticated = move on 
		        		.switchIfEmpty(Mono.defer(()->{
		        			//.defer -> RUN code ONMY if someone makes a request 
		        			//NOT WHEN THE APP RUNS !! 
		        			if(jwtService.isTokenValid(jwt, userEmail)) {
		        				//creatign a BADGE! of veritication :
		        				UsernamePasswordAuthenticationToken authToken = new UsernamePasswordAuthenticationToken(
		                                userEmail, null, new ArrayList<>() // FUTURE PROOF! can add ROLES!!! 
		                                //i'll brag about this for the PFA
		                                //TODO brag
		                            );
		        				return chain.filter(exchange)//go to next filter BUT NI BADGE
		                                .contextWrite(ReactiveSecurityContextHolder.withAuthentication(authToken));
		        						//= go to next filter AND ATTACH THE verifiied badge! 
		        					//.contextWire = attaching user to the dataStream itself NOT A THREAD
		        			}
		        			return chain.filter(exchange);
		        		}));
			}catch(Exception e) {
				return chain.filter(exchange);
				//still allows BUT not verified 
				//-> .authenticated routes will still block
				//ISNTEAD OF A FREAKING CRASH nice
			}
	        //flatMap -> reutrns ONLY if there's data
	        //meaning if secuirt context is empty -> flatMap is SKIPPED
	        //switchIfEmpty = the else block of it :) 
	        		 
		}
		
}


