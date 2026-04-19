import { Component } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { Router,RouterLink } from '@angular/router';

@Component({
  selector: 'app-home',
  standalone: true,
  imports: [CommonModule, FormsModule,RouterLink],
  templateUrl: './home.html',
  styleUrl: './home.css'
})
export class HomeComponent {

  hydrophobicityMin: number | null = null;
  hydrophobicityMax: number | null = null;
  massMin: number | null = null;
  massMax: number | null = null;
  stabilityMin: number | null = null;
  bindingMin: number | null = null;

  sequenceName: string = '';
  animationText: string = '';
  generated: boolean = false;

   constructor(private router: Router) {}

  generate(): void {
    this.generated = true;
    this.animationText = 'Génération en cours...';
    setTimeout(() => {
      this.animationText = '✓ Séquence générée avec succès !';
    }, 1500);
  }

  analyser(): void {
    if (!this.generated) {
      alert('Veuillez d\'abord générer une séquence.');
      return;
    }
    alert('Analyse lancée !');
  }

    compare(): void {
    this.router.navigate(['/compare']);
  }
}