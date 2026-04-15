import { Routes } from '@angular/router';
import { AuthComponent } from './auth/auth.component';
import { DashboardComponent } from './dashboard/dashboard.component';
import { GeneratorComponent } from './generator/generator.component';
import { ComparisonComponent } from './comparison/comparison.component';

export const routes: Routes = [
  { path: '', redirectTo: '/auth', pathMatch: 'full' },
  { path: 'auth', component: AuthComponent },
  { path: 'dashboard', component: DashboardComponent },
  { path: 'generator', component: GeneratorComponent },
  { path: 'comparison', component: ComparisonComponent }
];