import { Component } from '@angular/core';
import { FormsModule } from '@angular/forms';
import { RouterLink } from '@angular/router';
import { Router } from '@angular/router';

@Component({
  selector: 'app-login',
  standalone: true,
  imports: [FormsModule, RouterLink],
  templateUrl: './login.html',
  styleUrl: './login.css'
})
export class LoginComponent {

  email: string = '';
  password: string = '';

  constructor(private router: Router) {}

  onLogin(): void {
    if (this.email && this.password) {
      this.router.navigate(['/home']);
    } else {
      alert('Veuillez remplir tous les champs.');
    }
  }
}