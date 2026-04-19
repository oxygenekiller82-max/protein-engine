import { Component } from '@angular/core';
import { CommonModule } from '@angular/common';
import { RouterLink } from '@angular/router';
import { FormsModule } from '@angular/forms';
import { Router } from '@angular/router';

@Component({
  selector: 'app-signup',
  standalone: true,
  imports: [CommonModule, RouterLink, FormsModule],
  templateUrl: './signup.html',
  styleUrl: './signup.css'
})
export class SignupComponent {

  message: string = '';

  user = {
    firstName: '',
    lastName: '',
    email: '',
    password: ''
  };

  confirmPassword: string = '';

  constructor(private router: Router) {}

  register(): void {
    if (!this.user.firstName || !this.user.lastName || !this.user.email || !this.user.password) {
      this.message = 'Veuillez remplir tous les champs.';
      return;
    }
    if (this.user.password !== this.confirmPassword) {
      this.message = 'Passwords do not match.';
      return;
    }
    this.message = '';
    console.log(this.user);
    this.router.navigate(['/home']);
  }
}