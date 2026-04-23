package RealStuffs.AA_PFA.services;

import java.time.Instant;
import java.util.Date;
import java.util.HashMap;
import java.util.Map;
import java.util.function.Function;

import javax.crypto.SecretKey;

import org.springframework.beans.factory.annotation.Value;
import org.springframework.security.core.userdetails.UserDetails;
import org.springframework.stereotype.Service;
import org.springframework.transaction.annotation.Transactional;

import RealStuffs.AA_PFA.model.RefreshToken;
import RealStuffs.AA_PFA.model.User;
import RealStuffs.AA_PFA.repositories.RefreshTokenRepository;
import RealStuffs.AA_PFA.repositories.UserRepository;
import io.jsonwebtoken.Claims;
import io.jsonwebtoken.Jwts;
import io.jsonwebtoken.io.Decoders;
import io.jsonwebtoken.security.Keys;
import lombok.AllArgsConstructor;
import lombok.RequiredArgsConstructor;

import jakarta.annotation.PostConstruct;

//there is a new way to make tokens!! 
//Keys.hmacKeyFor() instead of strings huh..
@Service
@Transactional
@RequiredArgsConstructor
public class JwtService {
	
	private final RefreshTokenRepository refreshTokenRepository;
    private final UserRepository userRepository;
	
	//i forgor how to pull those from applciation.properties
	@Value("${application.security.jwt.secret-key}")
	private String secretKey;
	
	@Value("${application.security.jwt.expiration}")
    private long jwtExpiration;
	
	@Value("${application.security.jwt.refresh-token.expiration}")
    private long refreshExpiration;
	


	@PostConstruct
	public void debugProps() {
	    System.out.println("=== JWT CONFIG ===");
	    System.out.println("Expiration ms: " + jwtExpiration);
	    System.out.println("Refresh ms: " + refreshExpiration);
	    System.out.println("==================");
	}
	
	
	

	
	
	private SecretKey getSigningKey() {
		//striing key -> HMAC key -> bytes array 
        byte[]keyBytes= Decoders.BASE64.decode(secretKey);
        return Keys.hmacShaKeyFor(keyBytes);
    }
	
	//Generic extract claims (fields) from a toekn, jsut like servlet 
	public <T> T extractClaim(String token, Function<Claims, T> claimsResolver) {
	    final Claims claims = extractAllClaims(token);
	    return claimsResolver.apply(claims);
	}
	
	//extracting all claims for the token 
	private Claims extractAllClaims(String token) {
	    return Jwts.parser()
	    		.clockSkewSeconds(60) // WHAT
	            .verifyWith(getSigningKey()) // key from app properites btw
	            .build()
	            .parseSignedClaims(token)
	            .getPayload();
	}
	

	//token -> get userName 
	public String extractUsername(String token) {
	    return extractClaim(token, Claims::getSubject);
	}
	
	//token checker
	public boolean isTokenValid(String token, String userEmail) {
	    final String username = extractUsername(token);
	    //bug ?
	    
	    Date expiration = extractExpiration(token);
	    
	    return (username.equals(userEmail)) && !isTokenExpired(token);
	}
	
	private boolean isTokenExpired(String token) {
	    return extractExpiration(token).before(new Date());
	}

	private Date extractExpiration(String token) {
	    return extractClaim(token, Claims::getExpiration);
	}
	
	
	//finally genrating the stupid tokens 
	private String buildToken(
			Map<String, Object> extraClaims, 
			//FUTURE PROOFING Role adidng TODO brag
	        UserDetails userDetails, 
	        long expiration
			) {
		return Jwts.builder()
	            .claims(extraClaims)
	            .subject(userDetails.getUsername())
	            .issuedAt(new Date(System.currentTimeMillis()))
	            .expiration(new Date(System.currentTimeMillis() + expiration))
	            .signWith(getSigningKey()) // This uses your secret from properties
	            .compact();
	}
	
	
	//access token 
	public String generateToken(UserDetails userDetails) {
	    return buildToken(new HashMap<>(), userDetails, jwtExpiration);
	}
	
	//refresh 
	//THIS ALSO SAVES
	public String generateRefreshToken(UserDetails userDetails) {
	    String refreshToken= buildToken(new HashMap<>(), userDetails, refreshExpiration);
	    
	    //who ?
	    User user = userRepository.findByEmail(userDetails.getUsername())
	            .orElseThrow(() -> new RuntimeException("User not found"));
	    
	    //already has a token ?
	    RefreshToken tokenEntity = refreshTokenRepository.findByUser(user)
	    		.orElse(new RefreshToken());
	    
	    
	    RefreshToken refreshTokenEntity = new RefreshToken();
	    tokenEntity.setUser(user);
	    tokenEntity.setToken(refreshToken); 
	    tokenEntity.setExpiryDate(Instant.now().plusMillis(refreshExpiration));
	    
	    refreshTokenRepository.save(tokenEntity);
	    
	    return refreshToken;
	    
	}
	
	//stupid refresh 
	public String refreshAccessToken(String refreshToken) {
		//find current token 
		return refreshTokenRepository.findByToken(refreshToken)
	            .map(tokenEntity -> {
	                // Check if the Instant in the DB is before right now
	                if (tokenEntity.getExpiryDate().isBefore(Instant.now())) {
	                    refreshTokenRepository.delete(tokenEntity);
	                    throw new RuntimeException("Refresh token expired. Please log in again.");
	                }
	                return tokenEntity;
	                
	            })
	            
	            //if good new token 
	            .map(tokenEntity -> {
	                User user = tokenEntity.getUser();
	             
	                return generateToken(user); 
	            })
	            .orElseThrow(() -> new RuntimeException("Invalid Refresh Token. Not found in database."));
	            
       
	}
}
