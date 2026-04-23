import { Component, inject } from '@angular/core';
import { CommonModule } from '@angular/common';
import { Router } from '@angular/router';
import { RouterLink } from '@angular/router';
import { FormsModule } from '@angular/forms';

//our service
import { Auth } from '../services/auth';

@Component({
  selector: 'app-signup',
  standalone: true,
  imports: [CommonModule, RouterLink, FormsModule],
  templateUrl: './signup.html',
  styleUrl: './signup.css'
})

export class SignupComponent { 
  private authService = inject(Auth);
  private router = inject(Router);
  message: string = "";

  user = {
    firstName: '',
    lastName: '',
    email: '',
    password: ''
  };

confirmPassword: string = "";
register() {
  //validation client side oo

  if (this.user.password !== this.confirmPassword) {
    this.message = "Passwords do not match ";
    return;
  }

  this.message = "Connecting to server... ";

  //NOW THE CHANGES! -> springboot 
  this.authService.register(this.user).subscribe({
    next:(response)=>{
      console.log("Spring says:", response);
        this.message = "Account created, next page.";

        //delayt to see message
        setTimeout(() => {
          this.router.navigate(['/login']);
        }, 2000);
    },
    error:(err)=>{
      console.error("Spring rejected :", err);
      this.message = "Error: Email might already be taken or server is down";
    }
  });
}
}