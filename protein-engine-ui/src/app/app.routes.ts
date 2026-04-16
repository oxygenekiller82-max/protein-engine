import { Routes } from '@angular/router';
import { LoginComponent } from './login/login';      /* maktoub login 5ater ism el file login.ts */
import { SignupComponent } from './signup/signup';   /* maktoub signup 5ater ism el file signup.ts */

export const routes: Routes = [
  { path: 'login', component: LoginComponent },
  { path: 'signup', component: SignupComponent },
  { path: '', redirectTo: 'login', pathMatch: 'full' }
];