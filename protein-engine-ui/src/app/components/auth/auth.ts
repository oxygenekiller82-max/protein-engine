import { Component } from '@angular/core';
import { Router } from '@angular/router';

@Component({
  selector: 'app-auth',
  templateUrl: './auth.html',
  styleUrls: ['./auth.css']
})
export class AuthComponent {
  email = '';
  password = '';

  constructor(private router: Router) {}

  login() {
    if (this.email && this.password) {
      localStorage.setItem('token', 'fake-jwt-token');
      this.router.navigate(['/dashboard']);
    } else {
      alert('Entrez email et mot de passe');
    }
  }
}