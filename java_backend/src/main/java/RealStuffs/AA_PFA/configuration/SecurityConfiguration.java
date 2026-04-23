package RealStuffs.AA_PFA.configuration;

import java.util.List;

import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import org.springframework.http.HttpMethod;
import org.springframework.security.config.annotation.web.reactive.EnableWebFluxSecurity;
import org.springframework.security.config.web.server.SecurityWebFiltersOrder;
import org.springframework.security.config.web.server.ServerHttpSecurity;
import org.springframework.security.crypto.bcrypt.BCryptPasswordEncoder;
import org.springframework.security.web.server.SecurityWebFilterChain;
import org.springframework.security.web.server.context.NoOpServerSecurityContextRepository;
import org.springframework.web.cors.CorsConfiguration;
import org.springframework.web.cors.reactive.CorsConfigurationSource;
import org.springframework.web.server.ServerWebExchange;


//bean to inject for Bcrpyt 
@Configuration
@EnableWebFluxSecurity 
public class SecurityConfiguration {
	private final JwtAuthenticationFIlter jwtAuthFilter;
	
	@Bean
    public BCryptPasswordEncoder passwordEncoder() {
        return new BCryptPasswordEncoder();
    }
	
	public SecurityConfiguration(JwtAuthenticationFIlter jwtAuthFilter) {
        this.jwtAuthFilter = jwtAuthFilter;
    }
	
	@Bean 
	public SecurityWebFilterChain springSecurityFilterChain(ServerHttpSecurity http) {
		return http
				.csrf(ServerHttpSecurity.CsrfSpec::disable)
				.cors(cors -> cors.configurationSource(corsConfigurationSource()))
				.authorizeExchange(exchanges -> exchanges
						.pathMatchers(HttpMethod.OPTIONS, "/**").permitAll()
		                .pathMatchers("/api/auth/**","/api/auth/refresh").permitAll() // Login sign up for everyone
		                .pathMatchers("/api/**").permitAll()
		                .anyExchange().authenticated() // else give token hol up!
		            )
				
				.addFilterAt(jwtAuthFilter, SecurityWebFiltersOrder.AUTHENTICATION)
				//no more SESSIONS in reactive.. what
				.securityContextRepository(NoOpServerSecurityContextRepository.getInstance())
				//-> everyone is  stranger UNLESS they give a JWT+ no sessions
				.build(); // SECUTIY REPO IS CRAZY
				
	}
	
	//cors function 
	private CorsConfigurationSource corsConfigurationSource() {
        return (ServerWebExchange exchange) -> {
            CorsConfiguration config = new CorsConfiguration();
            config.setAllowedOrigins(List.of("http://localhost:4200"));//angular pls 
            config.setAllowedMethods(List.of("GET", "POST", "PUT", "DELETE", "OPTIONS"));
            config.setAllowedHeaders(List.of("*"));
            return config;
        };
    }
	
	

}
