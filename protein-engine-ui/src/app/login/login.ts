import { Component } from '@angular/core';
import { CommonModule } from '@angular/common';
import { RouterModule } from '@angular/router';

@Component({
  selector: 'app-login',
  standalone: true,
  imports: [CommonModule, RouterModule],
  templateUrl: './login.html',    /* Ism el file 3andek login.html */
  styleUrls: ['./login.css']       /* Ism el file 3andek login.css */
})
export class LoginComponent { }