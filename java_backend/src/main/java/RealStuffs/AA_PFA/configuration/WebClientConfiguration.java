package RealStuffs.AA_PFA.configuration;

import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import org.springframework.security.config.annotation.web.reactive.EnableWebFluxSecurity;
import org.springframework.security.config.web.server.ServerHttpSecurity;
import org.springframework.security.web.server.SecurityWebFilterChain;
import org.springframework.web.reactive.function.client.WebClient;

@Configuration
@EnableWebFluxSecurity
//to talk to python 
//@configuration -> when spring starts up it reads this file 
//BEFORE it handles any requests
public class WebClientConfiguration {
	@Bean 
	//bean objects -> can do @Autowired ! ->can call 
	//this WebCLinet in ANY SERVICE/controllers
	//WITHOUT doing new WebClient() !!
	public WebClient flaskWebClient() {
		return WebClient.builder()
				.baseUrl("http://localhost:5000")
				
				.codecs(configurer -> configurer
		                .defaultCodecs()
		                .maxInMemorySize(512 * 1024 * 1024)) //512MB bucket for that API..
				.build();
	}
	// spring already gives WebClient.Builder 
	//Builder = a tempalate -> u can customize it 
	//uri -> can just do .uri("/generate") ! -> goes to that uri
	//.build -> takes all settings(URL, headers, timeouts..)
	//-> freezes them in a single obeject -> WebClient ! 
	
	//POSTMAN ACCESS DENIED FIX, secuity blocks URL so..

	@Bean
	public SecurityWebFilterChain springSecurityFilterChain(ServerHttpSecurity http)  {
		    return http
		        .csrf(csrf -> csrf.disable()) 
		        .authorizeExchange(exchanges -> exchanges
		            .anyExchange().permitAll() //now good! 
		        )
		    	.build();
	}
	
}
