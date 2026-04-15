import { Component } from '@angular/core';
import { Router } from '@angular/router';
import { ApiService } from '../services/api.service';

@Component({
  selector: 'app-generator',
  templateUrl: './generator.html',
  styleUrls: ['./generator.css']
})
export class GeneratorComponent {
  target_length = 60;
  hydro_min = -15;
  hydro_max = 5;
  mass_min = 6200;
  mass_max = 6800;
  binding_min = 120;
  stability_min = 75;
  biological_switch = true;
  loading = false;
  result: any = null;

  constructor(private router: Router, private api: ApiService) {
    if (!localStorage.getItem('token')) {
      this.router.navigate(['/auth']);
    }
  }

  generate() {
    this.loading = true;
    const payload = {
      user_targets: {
        hydro_min: this.hydro_min,
        hydro_max: this.hydro_max,
        mass_min: this.mass_min,
        mass_max: this.mass_max,
        binding_min: this.binding_min,
        stability_min: this.stability_min
      },
      target_length: this.target_length,
      biological_switch: this.biological_switch
    };

    this.api.startSearch(payload).subscribe({
      next: (res) => {
        this.result = res;
        this.loading = false;
      },
      error: () => {
        this.result = { status: 'error', message: 'Erreur de connexion au serveur' };
        this.loading = false;
      }
    });
  }

  logout() {
    localStorage.removeItem('token');
    localStorage.removeItem('user');
    this.router.navigate(['/auth']);
  }
}