import { Routes } from '@angular/router';
import { LoginComponent } from './login/login';
import { SignupComponent } from './signup/signup';
import { HomeComponent } from './home/home';
import { CompareComponent } from './compare/compare';
import { LogsComponent } from './logs/logs';
import { AboutComponent } from './about/about';

export const routes: Routes = [
  { path: 'login', component: LoginComponent },
  { path: 'signup', component: SignupComponent },
  { path: 'home', component: HomeComponent },
  { path: 'compare', component: CompareComponent },
  { path: 'logs', component: LogsComponent },
  { path: 'about', component: AboutComponent },
  { path: '', redirectTo: 'login', pathMatch: 'full' }
];