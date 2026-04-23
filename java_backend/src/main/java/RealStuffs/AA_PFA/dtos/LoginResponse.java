package RealStuffs.AA_PFA.dtos;

import lombok.AllArgsConstructor;
import lombok.Builder;
import lombok.Data;
import lombok.NoArgsConstructor;

@Data
@Builder
@AllArgsConstructor
@NoArgsConstructor
public class LoginResponse {
	private String accessToken;
	private String refreshToken;
	
    private String username;
}
