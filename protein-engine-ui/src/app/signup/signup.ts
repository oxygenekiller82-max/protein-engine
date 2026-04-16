import { Component } from '@angular/core';
import { CommonModule } from '@angular/common';
import { RouterModule } from '@angular/router';
import { FormsModule } from '@angular/forms';
@Component({
  selector: 'app-signup',
  standalone: true,
  imports: [CommonModule, RouterModule,FormsModule],
  templateUrl: './signup.html',   /* Ism el file 3andek signup.html */
  styleUrl: './signup.css'       /* Ism el file 3andek signup.css */
})
export class SignupComponent { 
  message: string = "";

user = {
  firstName: '',
  lastName: '',
  email: '',
  password: ''
};

confirmPassword: string = "";
register() {

  if (this.user.password !== this.confirmPassword) {
    this.message = "Passwords do not match ";
    return;
  }

  this.message = "Account created successfully! ";

  console.log(this.user);
}
}
