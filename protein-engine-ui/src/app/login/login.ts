import { Component, inject,  } from '@angular/core';
import { FormsModule } from '@angular/forms';
import { RouterLink } from '@angular/router';
import { Router } from '@angular/router';
import { CommonModule } from '@angular/common'; // -> ngIf
import { Auth } from '../services/auth';

@Component({
  selector: 'app-login',
  standalone: true,
  imports: [FormsModule, RouterLink,CommonModule],
  templateUrl: './login.html',
  styleUrl: './login.css'
})
export class LoginComponent {
  //injection -> modern = inkec() oh
  private authService = inject(Auth);
  private router = inject(Router);


  email: string = '';
  password: string = '';

  message: string = '';


  onLogin(): void {
    if (!this.email || !this.password) {
      this.message=('Please do fill in all the fields');
      return;
      
    } 
      
    this.message = 'Connecting to server..';

    //Crendetials Object -> for auth
    const credentials = {
      email: this.email,
      password: this.password
    };

    //BACKEND CALL! 
    this.authService.login(credentials).subscribe({
      next: (response)=>{
        console.log('Success:', response);

        //JWT in localStorage
        localStorage.setItem('accessToken', response.accessToken);
        localStorage.setItem('refreshToken', response.refreshToken);
        localStorage.setItem('username', response.username);

        this.message = `Welcome back, ${response.username}!`;

        setTimeout(() => {
          this.router.navigate(['/home']);
        }, 1200);
      },
      error: (err) => {
        console.error('Login Error:', err);
        // spring login errors 
        this.message = 'Invalid email or password.';
      }
      });
  }
}